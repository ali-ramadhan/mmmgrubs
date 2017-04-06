library('FME')

OCS <- function (pars, constants) {
  derivs <- function(time, y, pars) {
    with (as.list(c(pars, constants, y)), {
      # Physical constants
      amu <- 1.66053886e-27 # [kg], atomic mass unit or 1 Dalton (Da)
      q0  <- 1.60217646e-19 # [C], elementary charge
      k   <- 8.987551e9     # [N m^2 / C^2], electrostatic constant
      
      m1 <- amu*m1Da; m2 <- amu*m2Da; m3 <- amu*m3Da
      q1 <- q0*q1e;   q2 <- q0*q2e;   q3 <- q0*q3e
      
      r12 <- sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
      r13 <- sqrt( (x1-x3)^2 + (y1-y3)^2 + (z1-z3)^2 )
      r23 <- sqrt( (x2-x3)^2 + (y2-y3)^2 + (z2-z3)^2 )
      
      # Position derivatives
      dx1 <- px1/m1;
      dy1 <- py1/m1;
      dz1 <- pz1/m1;
      
      dx2 <- px2/m2;
      dy2 <- py2/m2;
      dz2 <- pz2/m2;
      
      dx3 <- px3/m3;
      dy3 <- py3/m3;
      dz3 <- pz3/m3;
      
      # Momentum derivatives
      dpx1 <- k*q1*( q2*(x1-x2)/r12^3 + q3*(x1-x3)/r13^3 )
      dpy1 <- k*q1*( q2*(y1-y2)/r12^3 + q3*(y1-y3)/r13^3 )
      dpz1 <- k*q1*( q2*(z1-z2)/r12^3 + q3*(z1-z3)/r13^3 )
      
      dpx2 <- k*q2*( q1*(x2-x1)/r12^3 + q3*(x2-x3)/r23^3 )
      dpy2 <- k*q2*( q1*(y2-y1)/r12^3 + q3*(y2-y3)/r23^3 )
      dpz2 <- k*q2*( q1*(z2-z1)/r12^3 + q3*(z2-z3)/r23^3 )
      
      dpx3 <- k*q3*( q1*(x3-x1)/r13^3 + q2*(x3-x2)/r23^3 )
      dpy3 <- k*q3*( q1*(y3-y1)/r13^3 + q2*(y3-y2)/r23^3 )
      dpz3 <- k*q3*( q1*(z3-z1)/r13^3 + q2*(z3-z2)/r23^3 )
      
      d <- list(c(dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, dpx1, dpy1, dpz1,
                  dpx2, dpy2, dpz2, dpx3, dpy3, dpz3))
      return(d)
    })
  }
  
  # Now [pm]
  r12 <- with(as.list(pars), r12*1e-12)
  r23 <- with(as.list(pars), r23*1e-12)
  theta <- with(as.list(pars), theta)

  # Place the first atom to the left of central atom.
  x1_0 <- -r12
  y1_0 <- 0
  z1_0 <- 0
    
  # Place the central atom at the origin.
  x2_0 <- 0; y2_0 <- 0; z2_0 <- 0
  
  # Place the third atom to the right of the central taking into the account
  # the angle between the two bond lengths.
  x3_0 <- r23 * cos( (180 - theta)*(pi/180) )
  y3_0 <- r23 * sin( (180 - theta)*(pi/180) )
  z3_0 <- 0
  
  # Zero initial momentum
  px1_0 <- 0; py1_0 <- 0; pz1_0 <- 0
  px2_0 <- 0; py2_0 <- 0; pz2_0 <- 0
  px3_0 <- 0; py3_0 <- 0; pz3_0 <- 0
  
  y <- c(x1 = x1_0, y1 = y1_0, z1 = z1_0, x2 = x2_0, y2 = y2_0, z2 = z2_0,
         x3 = x3_0, y3 = y3_0, z3 = z3_0, px1 = px1_0, py1 = py1_0,
         pz1 = pz1_0, px2 = px2_0, py2 = py2_0, pz2 = pz2_0, px3 = px3_0,
         py3 = py3_0, pz3 = pz3_0)
  times <- c(seq(0, 1e-17, 1e-18), seq(1e-17, 100e-17, 1e-17), seq(1e-15, 100e-15, 1e-15))
  
  out <- ode(y = y, parms = pars, times = times, func = derivs)
  
  as.data.frame(out)
}

# For the cost I'll just use
# Make sure that measuredP = (p1x, p1y, p2x).
OCScost <- function (pars) {
  constants <- c(m1Da = 16, m2Da = 12, m3Da = 32, q1e = 2, q2e = 2, q3e = 2)
  measuredP <- c(4.5667e-22, -6.0556e-22, 1.9972e-22)
  
  out <- OCS(pars, constants)
  
  # Extract asymptotic momenta (ignore z components as they're all zero)
  px1 <- tail(out$px1,1); py1 <- tail(out$py1,1); # pz1 <- tail(out$pz1,1)
  px2 <- tail(out$px2,1); py2 <- tail(out$py2,1); # pz2 <- tail(out$pz2,1)
  px3 <- tail(out$px3,1); py3 <- tail(out$py3,1); # pz3 <- tail(out$pz3,1)
  
  p1 <- c(px1, py1); p2 <- c(px2, py2); p3 <- c(px3, py3)
  
  # Calculate the angle between the central atom and the +x-axis then rotate
  # the three momentum vectors back towards the origin by that much so that
  # the central's momentum is always along the +x-axis.
  
  phi <- atan2(py2, px2)
  R <- matrix( c(cos(-phi), -sin(-phi), sin(-phi), cos(-phi)), nrow=2, ncol=2, byrow=TRUE)
  
  p1 <- R %*% p1; p2 <- R %*% p2; p3 <- R %*% p3
  
  # We should have that p3 = - (p2 + p1) by conservation of momentum so
  # we're just going to use throw it out basically and only use the
  # independent momentum values (p1x, p1y, p2x) in our phase 1 model.
  
  # pinf <- c(0,0,0)
  # pinf[1] <- p1[1]; pinf[2] <- p1[2]; pinf[2] <- p2[1]

  asymptoticTime <- tail(out$time,1)  
  
  out <- out[-nrow(out),] # Delete last row
  
  # Create new row with our "momentum" and put it in.
  newrow <- c(asymptoticTime, rep(0,9), p1[1]*1e22, p1[2]*1e22, 0, p2[1]*1e22, 0, 0, 0, 0, 0)
  out <- rbind(out, newrow)

  DataPx1 <- c(time = asymptoticTime, px1 = measuredP[1]*1e22)
  DataPy1 <- c(time = asymptoticTime, py1 = measuredP[2]*1e22)
  DataPx2 <- c(time = asymptoticTime, px2 = measuredP[3]*1e22)
  
  cost <- modCost(model = out, obs = DataPx1, err = NULL)
  cost <- modCost(model = out, obs = DataPy1, err = NULL, cost = cost)
  cost <- modCost(model = out, obs = DataPx2, err = NULL, cost = cost)
  
  return(cost)
}

pars <- c(r12 = 115, r23 = 156, theta=175)
constants <- c(m1Da = 16, m2Da = 12, m3Da = 32, q1e = 2, q2e = 2, q3e = 2)
out <- OCS(pars, constants)

measuredP <- c(4.5667e-22, -6.0556e-22, 1.9972e-22)

# Model cost
cost <- OCScost(pars)
print(cost$model)

# Fit model to simulated data
# OCScost2 <- function(lpars)
#   OCScost(c(exp(lpars)))

# Pars <- pars[1:3] * 2
# Fit <- modFit(f = OCScost2, p = log(Pars), lower = log(c(50, 50, 50)), upper = log(c(1000, 1000, 360)))
# print(exp(coef(Fit)))

# Fit model to data
Fit <- modFit(f = OCScost, p = pars)

var0 <- Fit$var_ms_unweighted
# cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
MCMC <- modMCMC(f = OCScost, p = Fit$par, niter = 25000, jump = cov0,
                var0 = var0, wvar0 = 0.1, updatecov = 5, verbose=TRUE)
# MCMC$pars <- exp(MCMC$pars)