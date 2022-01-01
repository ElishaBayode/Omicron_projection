functions {
  //real switch_eta(real t, real t1, real eta, real nu, real xi) {
  //  return(eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu))));
  // }
 
 
  real[] sir(real t, //time
  real[] state,  //state
  real[] theta, //parameter
  real[] x_r, // data -real
  int[] x_i) { // data integer
  
  real S =  state[1];
  real Er = state[2];
  real Em = state[3];
  real Ir = state[4];
  real Im = state[5];
  real R =  state[6];
  real V =  state[7];
  real Erv = state[8];
  real Emv = state[9];
  real Irv = state[10];
  real Imv = state[11];
  real Rv =  state[12];
  real W =   state[13];
  real Erw = state[14];
  real Emw = state[15];
  real Irw = state[16];
  real Imw = state[17];
  real Rw =  state[18];
  real N  = x_r[1];
 // real beta_r = x_r[2];
 // real beta_m = x_r[3];
  real epsilon_r = x_r[2];
  real epsilon_m = x_r[3];
  real v = x_r[4]; //vax rate:
  real ve = x_r[5]; //:
  real sigma = x_r[6];
  real c = x_r[7];   // change to sigmoid function
  real mu = x_r[8];
  real w1 = x_r[9]; // waning R to S
  real w2 = x_r[10]; // waning Rv to V 
  real w3 = x_r[11]; // waning Rw to W 
  real gamma = x_r[12]; // recovery rate
  real vaxlevel=  x_r[13];
  real port_wane=  x_r[14];
  real past_infection = x_r[15];
  real b = x_r[16]; // boosters
  real nu = x_r[17];
  real beta_r = theta[1];
  real beta_m = theta[2]; 
  real wf = theta[3];
  real dydt[18];
  real X;
  int day;
  day = 1;
  while ((day + 1) < floor(t)) day = day + 1;
  

  //N = S+Er+Em+Ir+Im+R+V+Erv+Emv+Irv+Imv+Rv+W+Erw+Emw+Irw+Imw+Rw;//total population 
//lambda_r = c*beta_r*(Ir + Irv + Irw); #force of infection resident strain
//lambda_m = c*beta_m*(Im + Imv + Imw); #force of infection mutant strain
    
 dydt[1]  =  mu*N - (c*beta_r*(Ir + Irv + Irw)+c*beta_m*(Im + Imv + Imw))*S/N  + w1*R -(mu + nu*ve)*S;
 dydt[2]  = c*beta_r*(Ir + Irv + Irw)*S/N + wf* epsilon_r*c*beta_r*(Ir + Irv + Irw)*R/N - (sigma+mu)*Er; 
 dydt[3]  = c*beta_m*(Im + Imv + Imw)*S/N + wf* epsilon_m*c*beta_m*(Im + Imv + Imw)*R/N - (sigma+mu)*Em;
 dydt[4]  = sigma*Er - (gamma + mu)*Ir;
 dydt[5]  = sigma*Em - (gamma + mu)*Im;
 dydt[6]  = gamma*(Ir + Im) - wf*(epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*R/N - (mu + w1)*R;
 dydt[7]  = nu*ve*S + w2*Rv+w3*W - (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*V/N - (mu + b*ve)*V; 
 dydt[8]  = epsilon_r*c*beta_r*(Ir + Irv + Irw)*V/N + wf*epsilon_r*c*beta_r*(Ir + Irv + Irw)*Rv/N - (sigma+mu)*Erv; 
 dydt[9]  = epsilon_m*c*beta_m*(Im + Imv + Imw)*V/N +wf* epsilon_m*c*beta_m*(Im + Imv + Imw)*Rv/N - (sigma+mu)*Emv; 
 dydt[10]  =  sigma*Erv - (gamma + mu)*Irv;
 dydt[11]  = sigma*Emv - (gamma + mu)*Imv;
 dydt[12]  =  gamma*(Irv + Imv) -wf* (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*Rv/N - (mu + w2)*Rv;
 dydt[13]  =   b*ve*V + w2*Rw - (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*W/N -(mu+ w3)*W;
 dydt[14]  = epsilon_r*c*beta_r*(Ir + Irv + Irw)*W/N + wf*epsilon_r*c*beta_r*(Ir + Irv + Irw)*Rw/N - (sigma+mu)*Erw; 
 dydt[15]  = epsilon_m*c*beta_m*(Im + Imv + Imw)*W/N + wf*epsilon_m*c*beta_m*(Im + Imv + Imw)*Rw/N - (sigma+mu)*Emw; 
 dydt[16]  = sigma*Erw - (gamma + mu)*Irw;
 dydt[17]  = sigma*Emw - (gamma + mu)*Imw;
 dydt[18]  = gamma*(Irw + Imw) - wf*(epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*Rw/N - (mu + w2)*Rw;
 
 return dydt;
  }
}

data{
int<lower=1> T;     // number of time steps
  int<lower=1> N;     // number of days
  int<lower=1> J;     // number of response data timeseries
  real y0_vars[17];   // initial state  
  real t0;            // first time step
  real time[T];       // time increments
  int days[N];        // day increments
  int last_day_obs;   // last day of observed data; days after this are projections 
  int daily_cases_res[last_day_obs,J]; // daily new case counts (resident)
  int daily_cases_mut[last_day_obs,J]; // daily new case counts (mutant)
  int daily_cases_tot[last_day_obs,J]; // daily new cases (total)
  int n_x_r;         // the number of x_r values
  real x_r[n_x_r];   // data for ODEs (real numbers)
  int n_x_i;         // the number of x_i values
  int x_i[n_x_i];
  real delay_scale[J];    // Weibull parameter for delay in becoming a case
  real delay_shape[J];    // Weibull parameter for delay in becoming a case
  int time_day_id[N]; // last time increment associated with each day
  //int time_day_id0[N];// first time increment for Weibull integration of cases
  real beta_r_prior[2];   // lognormal log mean and SD for R0 prior
  real i0r_prior[2]; 
  real beta_m_prior[2];   // lognormal log mean and SD for R0 prior
  real i0m_prior[2]; 
  real samp_frac_prior[2];  
  real phi_prior;     // SD of normal prior on 1/sqrt(phi) [NB2(mu, phi)]
  int use_phi_prior2;  // Logical
  real phi_prior2[2];   // alternate lognormal
  int<lower=0, upper=1> priors_only;
  int<lower=0, upper=J> est_phi; // estimate NB phi?
  int<lower=0, upper=1> obs_model; // observation model: 0 = Poisson, 1 = NB2
  real<lower=0> rw_sigma; // specified random walk standard deviation
  int<lower=0, upper=1> contains_NAs; // Logical: contains NA values?
  int samp_frac_seg[N]; // optional index of estimated sample fractions for 1st timeseries
  int<lower=1, upper=4> samp_frac_type;
  real samp_frac_fixed[N,J];  
  real ode_control[3]; // vector of ODE control numbers
  int<lower=0> K;      // number of linear predictors
  matrix[N,K] X; 
}

parameters {
  real<lower=0, upper=x_r[1]> i0r; 
  real<lower=0, upper=x_r[1]> i0m; 
  real beta_r; // Stan ODE solver seems to be more efficient without this bounded at > 0
  real beta_m; 
 real<lower=0> phi[est_phi]; // NB2 (inverse) dispersion; `est_phi` turns on/off
 real beta[K]; // 
  
}

transformed parameters {
  real dx = time[2] - time[1]; // time increment
  real ftr[T]; // container for the lambda function at time t
  real ftm[T];
  real ftot[T];
  real mur[N,J]; // estimated daily cases for each day
  real mum[N,J];
  real mutot[N,J];
  real sum_ftr_inner; // holds a temporary calculation
  real sum_ftm_inner;
  real sum_ftot_inner;
  real etar[N,J]; // expected value on link scale (log)
  real etam[N,J];
  real etatot[N,J];
  real sigma;
  real ve;
  real Er; 
  real Em; 
  real Erv; 
  real Erw; 
  real Emv; 
  real Emw; 
  real y_hat[T,18];  
  real this_samp; 
  real N_pop = x_r[1];
  real vaxlevel;
  real port_wane;
  real past_infection;
  real c;
  real incres; //incidence of resident strain on day 1
  real incmut; //incidence of mutant strain on day 1
  real epsilon_r;
  real epsilon_m;
  real ff; 
  real Vtot;
  real Wtot;
  real Ertot;
  real Emtot;
  real Ervw;
  real Emvw;
  //real this_samp; // holds the sample fraction for a given day
  real y0[18]; //
  Vtot = vaxlevel*N_pop*(1-port_wane);//allocate to V, Ev, Iv
  Wtot = vaxlevel*N*port_wane; //allocate to W, Ew, Iw 
  Ertot = ff*incres/sigma; 
  Emtot = ff*incmut/sigma;
  Ervw =  vaxlevel*ve*Ertot;
  Emvw =  vaxlevel*ve*Emtot;
 y0[1] = N_pop - (y0_vars[1]+y0_vars[2]+ y0_vars[3]+ y0_vars[4]+y0_vars[5]+y0_vars[6] + y0_vars[7]+y0_vars[8]+y0_vars[9]+y0_vars[10]+y0_vars[11]+y0_vars[12] + y0_vars[13]+y0_vars[14]+y0_vars[15]+y0_vars[16]+y0_vars[17]);
 y0[2] = (1-vaxlevel*epsilon_r)*Ertot;
 y0[3] = (1-vaxlevel*epsilon_m)*Emtot;
 y0[4] = ff*2*y0_vars[2];
 y0[5] = ff*2*y0_vars[3]; 
 y0[6] = N_pop*past_infection*(1-vaxlevel);
 y0[7] = N_pop *(1-past_infection)*vaxlevel*(1-port_wane) - (y0_vars[8]+y0_vars[10]+y0_vars[9]+y0_vars[11])*(1-port_wane);
 y0[8] = Ervw *(1-port_wane);
 y0[9] = Emvw *(1-port_wane);
 y0[10] = ff*2*y0_vars[8];
 y0[11] = ff*2*y0_vars[9];
 y0[12] = N_pop*past_infection*(vaxlevel) * (1-port_wane);
 y0[13] = N_pop*(1-past_infection)*vaxlevel*(port_wane) - (y0_vars[14]+y0_vars[15]+y0_vars[16]+y0_vars[17])*(port_wane); 
 y0[14] = Ervw * port_wane;
 y0[15] = Emvw * port_wane;
 y0[16] = ff*2*y0_vars[14];
 y0[17] = ff*2*y0_vars[15]; 
 y0[18] = N_pop*past_infection*(vaxlevel) * port_wane;
 real theta
 theta[1] = beta_r;
 theta[2] = beta_m;
 theta[3] = wf;

 //R0r = theta[1];
 //R0m = theta[2];
 
  y_hat = integrate_ode_rk45(sir, y0, t0, time, theta, x_r, x_i, ode_control[1], ode_control[2], ode_control[3]);


for (j in 1:J) { // data_type increment
    for (n in 1:N) { // day increment
      this_samp = samp_frac_fixed[n,j];
      if (n_samp_frac > 1 && j == 1) {
        if (samp_frac_type != 4) { // anything but segmented samp_frac
          if (n <= last_day_obs) {
            this_samp = samp_frac[n];
          }
          if (n > last_day_obs && j == 1) {
            this_samp = samp_frac[n_samp_frac]; // forecast with last value
          }
        } else { // segmented
          this_samp = samp_frac[samp_frac_seg[n]];
        }
      }
      if (n_samp_frac == 1 && j == 1) {
        this_samp = samp_frac[1];
      }
      for (t in 1:T) {
        ftr[t] = 0; 
        //ftm[t] = 0;
       //ftot[t] = 0;// initialize across the full 1:T
      }
      // a fancy way of moving across a window of time:
      for (t in time_day_id0[n]:time_day_id[n]) { // t is an increment here
       
    
        sigma = x_r[8];
        Er = y_hat[t,2];
        Erv = y_hat[t,8];
        Erw = y_hat[t,14]; 

        ftr[t] = this_samp * sigma * (Er + Erv + Erw) *
        exp(weibull_lpdf(time[time_day_id[n]] - time[t] | delay_shape[j], delay_scale[j]));
      }
      
      for (t in 1:T) {
        ftm[t] = 0; 
        //ftm[t] = 0;
       //ftot[t] = 0;// initialize across the full 1:T
      }
      
      for (t in time_day_id0[n]:time_day_id[n]) { // t is an increment here
       
        sigma = x_r[8];
        Em = y_hat[t,3];
        Emv = y_hat[t,9];
        Emw = y_hat[t,15]; 

        ftm[t] = this_samp * sigma * (Em + Emv + Emw) *
        exp(weibull_lpdf(time[time_day_id[n]] - time[t] | delay_shape[j], delay_scale[j]));
      }
      
      for (t in 1:T) {
        ftot[t] = 0; 
        //ftm[t] = 0;
       //ftot[t] = 0;// initialize across the full 1:T
      }
      
      for (t in time_day_id0[n]:time_day_id[n]) { // t is an increment here
       
        sigma = x_r[8];
        Em = y_hat[t,3];
        Emv = y_hat[t,9];
        Emw = y_hat[t,15]; 
        Em = y_hat[t,3];
        Emv = y_hat[t,9];
        Emw = y_hat[t,15];

        ftot[t] = this_samp * sigma * (Em + Emv + Emw + Er + Erv + Erw) *
        exp(weibull_lpdf(time[time_day_id[n]] - time[t] | delay_shape[j], delay_scale[j]));
      }
      
      
      sum_ftr_inner = 0; // initialize
      sum_ftm_inner = 0;
      sum_ftot_inner = 0;
      
      for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
        sum_ftr_inner += ftr[t];
      }
      
      for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
        sum_ftm_inner += ftm[t];
      }
      
      for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
        sum_ftot_inner += ftot[t];
      }
      // trapezoid integration:
      mur[n,j] = 0.5 * dx *
      (ftr[time_day_id0[n]] + 2 * sum_ftr_inner + ftr[time_day_id[n]]);
      etar[n,j] = log(mur[n,j]);

      if (K > 0) { // has a linear predictor
        for (k in 1:K) etar[n,j] = etar[n,j] + X[n,k] * beta[k];
        mur[n,j] = exp(etar[n,j]);
      }
    
    
     mum[n,j] = 0.5 * dx *
      (ftm[time_day_id0[n]] + 2 * sum_ftm_inner + ftm[time_day_id[n]]);
      etam[n,j] = log(mum[n,j]);

      if (K > 0) { // has a linear predictor
        for (k in 1:K) etam[n,j] = etam[n,j] + X[n,k] * beta[k];
        mum[n,j] = exp(etam[n,j]);
      }
      
       mutot[n,j] = 0.5 * dx *
      (ftot[time_day_id0[n]] + 2 * sum_ftot_inner + ftot[time_day_id[n]]);
      etatot[n,j] = log(mutot[n,j]);

      if (K > 0) { // has a linear predictor
        for (k in 1:K) etatot[n,j] = etatot[n,j] + X[n,k] * beta[k];
        mutot[n,j] = exp(etatot[n,j]);
      }
    
     }
}
}



model {
  // priors:
  if (est_phi > 0 && obs_model == 1) { // NB2
    // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    // log(abs(D(expression(1/sqrt(x)), "x"))); log(0.5 * x^-0.5/sqrt(x)^2
    // log absolute derivative of the transform

    if (use_phi_prior2) {
      phi ~ lognormal(phi_prior2[1], phi_prior2[2]);
    } else {
      for (j in 1:J) {
        1/sqrt(phi[j]) ~ normal(0, phi_prior);
        target += log(0.5) - 1.5 * log(phi[j]); // Jacobian adjustment
      }
    }
  }
  beta ~ normal(0, 1);
  beta_r ~ lognormal(beta_r_prior[1], beta_r_prior[2]);
  beta_m ~ lognormal(beta_m_prior[1], beta_m_prior[2]);
  i0r ~ lognormal(i0_prior[1], i0_prior[2]);
  i0m ~ lognormal(i0_prior[1], i0_prior[2]);
  //fsi ~ beta(e_prior[1], e_prior[2]); // two names
  // log(abs(D(expression(ud / (ud + ur)), "ur"))); log(abs(-(ud/(ud + ur)^2))
  // log absolute derivative of the transform
 // target += log(ud) - 2 * log(ud + ur); // Jacobian adjustment

 // start_decline ~ lognormal(start_decline_prior[1], start_decline_prior[2]);
//  end_decline ~ lognormal(end_decline_prior[1], end_decline_prior[2]);
 // for (s in 1:S) {
 //   f_s[s] ~ beta(f_prior[s,1], f_prior[s,2]); // allow separate f priors
//  }
  if (n_samp_frac > 0 && samp_frac_type != 4) { // samp_frac estimated but not segmented
    samp_frac[1] ~ beta(samp_frac_prior[1], samp_frac_prior[2]);
    if (n_samp_frac > 1) {
      for (n in 2:n_samp_frac) {
        samp_frac[n] ~ normal(samp_frac[n - 1], rw_sigma); // RW
      }
    }
  }
  
  if (n_samp_frac > 0 && samp_frac_type != 4) { // samp_frac segmented
    for (n in 1:n_samp_frac) {
      samp_frac[n] ~ beta(samp_frac_prior[1], samp_frac_prior[2]);
    }
  }

  // data likelihood resident:
  
  
  if (!priors_only) { // useful to turn off for prior predictive checks
    if (contains_NAs) { // Not vectorized to 'easily' deal with NAs:
      for (n in 1:last_day_obs) {
        for (j in 1:J) {
          if (daily_cases_res[n,j] != 9999999) { // NA magic number
            if (obs_model == 0) {
              daily_cases_res[n,j] ~ poisson_log(etar[n,j]);
            } else if (obs_model == 1) {
              daily_cases_res[n,j] ~ neg_binomial_2_log(etar[n,j], phi[j]);
            }
          }
        }
      }
    } else { // No NAs; vectorized for increased efficiency:
      for (j in 1:J) {
        if (obs_model == 0) {
          daily_cases_res[1:last_day_obs,j] ~ poisson_log(etar[1:last_day_obs,j]);
        } else if (obs_model == 1) {
          daily_cases_res[1:last_day_obs,j] ~ neg_binomial_2_log(etar[1:last_day_obs,j], phi[j]);
        }
      }
    }
  }

if (!priors_only) { // useful to turn off for prior predictive checks
    if (contains_NAs) { // Not vectorized to 'easily' deal with NAs:
      for (n in 1:last_day_obs) {
        for (j in 1:J) {
          if (daily_cases_mut[n,j] != 9999999) { // NA magic number
            if (obs_model == 0) {
              daily_cases_mut[n,j] ~ poisson_log(etam[n,j]);
            } else if (obs_model == 1) {
              daily_cases_mut[n,j] ~ neg_binomial_2_log(etam[n,j], phi[j]);
            }
          }
        }
      }
    } else { // No NAs; vectorized for increased efficiency:
      for (j in 1:J) {
        if (obs_model == 0) {
          daily_cases_mut[1:last_day_obs,j] ~ poisson_log(etam[1:last_day_obs,j]);
        } else if (obs_model == 1) {
          daily_cases_mut[1:last_day_obs,j] ~ neg_binomial_2_log(etam[1:last_day_obs,j], phi[j]);
        }
      }
    }
  }
  
  
  if (!priors_only) { // useful to turn off for prior predictive checks
    if (contains_NAs) { // Not vectorized to 'easily' deal with NAs:
      for (n in 1:last_day_obs) {
        for (j in 1:J) {
          if (daily_cases_tot[n,j] != 9999999) { // NA magic number
            if (obs_model == 0) {
              daily_cases_tot[n,j] ~ poisson_log(etatot[n,j]);
            } else if (obs_model == 1) {
              daily_cases_tot[n,j] ~ neg_binomial_2_log(etatot[n,j], phi[j]);
            }
          }
        }
      }
    } else { // No NAs; vectorized for increased efficiency:
      for (j in 1:J) {
        if (obs_model == 0) {
          daily_cases_tot[1:last_day_obs,j] ~ poisson_log(etatot[1:last_day_obs,j]);
        } else if (obs_model == 1) {
          daily_cases_tot[1:last_day_obs,j] ~ neg_binomial_2_log(etatot[1:last_day_obs,j], phi[j]);
        }
      }
    }
  }
  
  
}
generated quantities{
  
  // real e; // renamed fsi
  int y_rep_res[N,J]; // posterior predictive replicates
  for (j in 1:J) {
    for (n in 1:N) {
      if (obs_model == 0) {
       y_rep_res[n,j] = poisson_log_rng(etar[n,j]);
      } else if (obs_model == 1) {
        y_rep_res[n,j] = neg_binomial_2_log_rng(etar[n,j], phi[1]);
      }
    }
  }
   int y_rep_mut[N,J]; // posterior predictive replicates
  for (j in 1:J) {
    for (n in 1:N) {
      if (obs_model == 0) {
        y_rep_mut[n,j] = poisson_log_rng(etam[n,j]);
      } else if (obs_model == 1) {
        y_rep_mut[n,j] = neg_binomial_2_log_rng(etam[n,j], phi[1]);
      }
    }
  }
  int y_rep_tot[N,J]; // posterior predictive replicates
  for (j in 1:J) {
    for (n in 1:N) {
      if (obs_model == 0) {
        y_rep_tot[n,j] = poisson_log_rng(etatot[n,j]);
      } else if (obs_model == 1) {
        y_rep_tot[n,j] = neg_binomial_2_log_rng(etatot[n,j], phi[1]);
      }
    }
  }
 real R0r = beta_r / gamma;
 real R0m = beta_m / gamma;
 
}  
