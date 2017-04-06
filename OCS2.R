library('rstan')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n <- 1
p <- c(4.5667e-22, -6.0556e-22, 1.9972e-22)

# times <- scan('times.txt')
# times <- as.list(read.table('times.txt'))


r12 <- 115.01
r23 <- 156.0
theta <- 170.0


x1 <- -r12*1e-12
x3 <- r23*1e-12 * cos( (pi/180)*(180 - theta) );
y3 <- r23*1e-12 * sin( (pi/180)*(180 - theta) );

dat <- list(x_r = c(16.01, 12.02, 32.03),
           x_i = c(2, 2, 2),
           t0 = 0,
           times = c(1e-19, 2e-19, 3e-19, 4e-19, 5e-19, 6e-19, 7e-19, 8e-19, 9e-19, 10e-19),
           pars = c(r12, r23, theta),
           q0 = c(x1, 0, 0, 0, 0, 0, x3, y3, 0, rep(0,9)) )

OCS_fit <- stan(file = 'OCS2.stan', sample_file="sample.log", data = dat,
                 diagnostic_file="diag.log", save_dso = TRUE, verbose = TRUE,
                 chains=1, iter = 10, warmup = 5)
