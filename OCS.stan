/* Stan model for infering the initial geometry of OCS in a Coulomb explosion
 * imaging experiment.
 */

functions {
  // Must always use this very specific signature, even if you don't use all
  // the arguments.
  real[] HamiltonsEqs(real t, real[] q, real[] theta, real[] x_r, int[] x_i) {
    /* Notes: q[1:9] = [x1, y1, z1, ..., x3, y3, z3]
     *        q[10:18] = [px1, py1, pz1, ..., px3, py3, pz3]
     */
    
    real dqdt[18]; 
    
    amu = 1.66053886e-27; // [kg], atomic mass unit
    e   = 1.60217646e-19; // [C], elementary charge
    k   = 8.987551e9;     // [N m^2 / C^2], electrostatic constant
    
    m1 = amu*x_r[1]; m2 = amu*x_r[2]; m3 = amu*x_r[3];
    q1 = e*x_i[1];   q2 = e*x_i[2];   q3 = e*x_i[3];
    
    r12 = sqrt( square(q[1]-q[4])  + square(q[2]-q[5])  + square(q[3]-q[6]) );
    r13 = sqrt( square(q[1]-q[7])  + square(q[2]-q[8])  + square(q[3]-q[9]) );
    r23 = sqrt( square(q[4]-q[7])  + square(q[5]-q[8])  + square(q[6]-q[9]) );
    
    // position derivatives
    dqdt[1] = q[10]/m1;
    dqdt[2] = q[11]/m1;
    dqdt[3] = q[12]/m1;
    
    dqdt[4] = q[13]/m2;
    dqdt[5] = q[14]/m2;
    dqdt[6] = q[15]/m2;
    
    dqdt[7] = q[16]/m3;
    dqdt[8] = q[17]/m3;
    dqdt[9] = q[18]/m3;
    
    // momentum derivatives
    dqdt[10] =  k*q1*( q2*(q[1]-q[4])/pow(r12,3) + q3*(q[1]-q[7])/pow(r13,3) );
    dqdt[11] =  k*q1*( q2*(q[2]-q[5])/pow(r12,3) + q3*(q[2]-q[8])/pow(r13,3) );
    dqdt[12] =  k*q1*( q2*(q[3]-q[6])/pow(r12,3) + q3*(q[3]-q[9])/pow(r13,3) );
    
    dqdt[13] =  k*q2*( q1*(q[4]-q[1])/pow(r12,3) + q3*(q[4]-q[7])/pow(r23,3) );
    dqdt[14] =  k*q2*( q1*(q[5]-q[2])/pow(r12,3) + q3*(q[5]-q[8])/pow(r23,3) );
    dqdt[15] =  k*q2*( q1*(q[6]-q[3])/pow(r12,3) + q3*(q[6]-q[9])/pow(r23,3) );
    
    dqdt[16] =  k*q3*( q1*(q[7]-q[1])/pow(r13,3) + q2*(q[7]-q[4])/pow(r23,3) );
    dqdt[17] =  k*q3*( q1*(q[8]-q[2])/pow(r13,3) + q2*(q[8]-q[5])/pow(r23,3) );
    dqdt[18] =  k*q3*( q1*(q[9]-q[3])/pow(r13,3) + q2*(q[9]-q[6])/pow(r23,3) );
    
    return dqdt;
  }
  
  real deg2rad(real deg) {
    return deg * (pi/180);
  }
  
  // ...
  real[3] coulombExplode(real r12, real r23, real theta) {
    // Place the first atom to the left of central atom.
    x1 = -r12;
    y1 = 0;
    z1 = 0;
    
    //Place the central atom at the origin.
    x2 = 0;
    y2 = 0;
    z2 = 0;
    
    // Place the third atom to the right of the central taking into the account
    // the angle between the two bond lengths.
    x3 = r23 * cos(deg2rad(180 - theta));
    y3 = r23 * sin(deg2rad(180 - theta));
    z3 = 0;
    
    g = [x_1 y_1 z_1 x_2 y_2 z_2 x_3 y_3 z_3];
    p_0 = zeros(1,9);
    
    y_hat = integrate_ode_bdf(sho, y0, t0, ts, theta, x_r, x_i,rel_tol, abs_tol, max_steps);
    
  }
}

data {
  int<lower=1> n; // number of CEI events observed (n=1 for now lol)
  real p[3]; // (p1x, p2x, p2y) measurements
}

transformed data {
  // Required empty values for the ODE solver function call.
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> r12; // [pm]
  real<lower=0> r23; // [pm]
  real<lower=0> theta; // [deg]
  /* theta should really have upper=180 but might introduce complications if I
   * use a normal prior for theta and end up samlping theta > 180. will come
   * back to this. Might not be a problem though as other models seem to do
   * something similar.
   */
   
   real<lower=0> sigma;
}

model {
  // "Priors" on latent variables from quantum chemistry simulations
  r12 ~ normal(114, 12.4);
  r23 ~ normal(163, 12.4);
  theta ~ normal(174.4, 3.4);
  
  // Weakly informative priors (give it a try)
  
  // Likelihood
  p_hat = CoulombExplode(r12, r23, theta);
  for (i in 1:n)
    p[i] ~ normal(p_hat, sigma);
  
  real y_hat[T,2];
  sigma ~ cauchy(0, 2.5);
  theta ~ normal(0, 1);
  y0 ~ normal(0, 1);
  y_hat = integrate_ode_rk45(sho, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T)
    y[t] ~ normal(y_hat[t], sigma);
}