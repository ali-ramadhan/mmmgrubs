library('rstan')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n <- 1
p <- c(4.5667e-22, -6.0556e-22, 1.9972e-22)
x_r <- c(16.0, 12.0, 32.0)
x_i <- c(2, 2, 2)
t0 <- 0
# times <- scan('times.txt')
# times <- as.list(read.table('times.txt'))
times <- c(1e-19, 2e-19, 3e-19, 4e-19, 5e-19, 6e-19, 7e-19, 8e-19, 9e-19, 10e-19)
r12 <- 115
r23 <- 156
theta <- 170

OCS_fit <- stan(file = 'OCS2.stan', sample_file="sample.log",
                 diagnostic_file="diag.log", save_dso = TRUE, verbose = TRUE,
                 chains=1, iter = 10, warmup = 5)
