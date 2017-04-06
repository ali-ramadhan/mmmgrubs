/* Stan model for infering the initial geometry of OCS in a Coulomb explosion
 * imaging experiment.
 */

functions {
  // Must always use this very specific signature, even if you don't use all
  // the arguments.
  real[] hamiltonsEqs(real t, real[] q, real[] theta, real[] x_r, int[] x_i) {
    /* Notes: q[1:9] = [x1, y1, z1, ..., x3, y3, z3]
     *        q[10:18] = [px1, py1, pz1, ..., px3, py3, pz3]
     */
    
    real dqdt[18];
    real amu; real e; real k;
    real m1; real m2; real m3;
    real q1; real q2; real q3;
    real r12; real r13; real r23;
    
    print("Entering hamiltonEqs()...");
    
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
    
    print("Leaving hamiltonEqs()...");
    return dqdt;
  }
  
  real deg2rad(real deg) {
    return deg * (pi()/180);
  }
  
  real[] coulombExplode(real r12, real r23, real theta,
                        real t0, real[] times, real[] x_r, int[] x_i) {
    /* Notes: q[1:9] = [x1, y1, z1, ..., x3, y3, z3]
     *        q[10:18] = [px1, py1, pz1, ..., px3, py3, pz3]
     */
    
    real q0[18];
    real q[10,18];
    real par[1];

/*    real rel_tol;
    real abs_tol;
    real max_num_steps;*/
    
    vector[3] p1; vector[3] p2; vector[3] p3;
    matrix[2,2] R;
    real theta_2x;
    vector[3] p;
    
    real r12m;
    real r23m;

    par[1] = 0;
    
    r12m = r12*1e-12; // [pm] -> [m]
    r23m = r23*1e-12; // [pm] -> [m]
    
    print("Entering coulombExplode()...");
    
    // Place the first atom to the left of central atom.
    q0[1] = -r12m;
    q0[2] = 0;
    q0[3] = 0;
    
    //Place the central atom at the origin.
    q0[4] = 0;
    q0[5] = 0;
    q0[6] = 0;
    
    // Place the third atom to the right of the central taking into the account
    // the angle between the two bond lengths.
    q0[7] = r23m * cos(deg2rad(180 - theta));
    q0[8] = r23m * sin(deg2rad(180 - theta));
    q0[9] = 0;
    
    // zero initial momentum
    q0[10] = 0; q0[11] = 0; q0[12] = 0;
    q0[13] = 0; q0[14] = 0; q0[15] = 0;
    q0[16] = 0; q0[17] = 0; q0[18] = 0;
    
/*    rel_tol = 1e-6;
    abs_tol = 1e-27;
    max_num_steps = 1000;*/
    
    q = integrate_ode_rk45(hamiltonsEqs, q0, t0, times, par, x_r, x_i);
/*    q = integrate_ode_rk45(hamiltonsEqs, q0, t0, times, par, x_r, x_i,
                              rel_tol, abs_tol, max_num_steps);*/
    
    // Extract each atom's momentum into 2D vectors in preparation to rotate.
    p1 = to_vector(q[5, 10:11]);
    p2 = to_vector(q[5, 13:14]);
    p3 = to_vector(q[5, 16:17]);
    
    // Put each momentum into column vector form so we can use vanilla matrix
    // multiplication.
    // p1 = p1'; p2 = p2'; p3 = p3';
    
    // Calculate the angle between the central atom and the +x-axis then rotate
    // the three momentum vectors back towards the origin by that much so that
    // the central's momentum is always along the +x-axis.
    theta_2x = atan2(p2[2], p2[1]);
    R[1,1] = cos(-theta_2x); R[1,2] = -sin(-theta_2x);
    R[2,1] = sin(-theta_2x); R[2,2] =  cos(-theta_2x);
    p1 = R*p1;
    p2 = R*p2;
    p3 = R*p3;
    
    // Put everything back into a row vector.
    // p1 = p1'; p2 = p2'; p3 = p3';
    
    // Set the z components (and also y in case of carbon) to zero so they all
    // have exactly the same value rather than 0.0000 and -0.0000, etc.
    // p1[3] = 0;
    // p2[2] = 0; p_2[3] = 0;
    // p3[3] = 0;
    
    // We should have that p3 = - (p2 + p1) by conservation of momentum so
    // we're just going to use throw it out basically and only use the
    // independent momentum values (p1x, p1y, p2x) in our phase 1 model.
    p[1] = p1[1]; p[2] = p1[2]; p[3] = p2[1];
    print("Leaving coulombExplode()...")
    
    return to_array_1d(p);
  }
}

data {
  int<lower=1> n; // number of CEI events observed (n=1 for now lol)
  real p[3]; // (p1x, p2x, p2y) measurements
  real<lower=0> x_r[3]; // masses of the three atoms in [amu].
  int<lower=1> x_i[3];  // charges of the three atoms in [e].
  real t0;
  real times[10]; // all the time steps for the ODE (I know...)
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

/*model {
  real p_hat[3];
  
  // "Priors" on latent variables from quantum chemistry simulations
  r12 ~ normal(114, 12.4);
  r23 ~ normal(163, 12.4);
  theta ~ normal(174.4, 3.4);
  
  // Weakly informative priors (give it a try?)
  
  // Likelihood
  p_hat = coulombExplode(r12, r23, theta, t0, times, x_r, x_i);
  p ~ normal(p_hat, sigma);
}*/

model {
}

generated quantities {
  real y_hat[3];
  
  print("Entering generated quantities...")
  y_hat <- coulombExplode(r12, r23, theta, t0, times, x_r, x_i);

/*  for (t in 1:T) {
    y_hat[t,1] <- y_hat[t,1] + normal_rng(0,sigma[1]);
    y_hat[t,2] <- y_hat[t,2] + normal_rng(0,sigma[2]);
  }*/
}
