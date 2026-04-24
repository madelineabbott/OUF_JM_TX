// fit joint model with tx effect modeled as drift and hazard form 2
functions {
  vector r_style_subset(matrix x, int col_id, array[] int row_ids) {
    vector[size(row_ids)] y;
    vector[rows(x)] x_col;
    int pos;
    pos = 1;
    x_col = col(x, col_id);
    for (i in 1:size(row_ids)) {
      y[pos] = x_col[row_ids[i]];
      pos = pos + 1;
    }
    return y;
  }
}

data {
	int Nall;					// Sample size (meas occ + grid)
	int Nlong;        // Sample size (meas occ only)
	int N;            // Number of subjects
	int Nevents;      // Number of recurrent events (obs. + censored)
  int Nall_tx;      // Sample size ((meas occ + grid) x num prior tx)	// 
	int K;					// Number of items
	int P;					// Number of latent factors

	array[Nlong] int meas_occ_rows; // for longitudinal measurement occasions
	array[N] int cumu_meas; // for longitudinal measurement occasions
	array[N] int repme_meas; // for longitudinal measurement occasions

	array[N] int cumu; // for longitudinal + grid measurement occasions
	array[N] int repme; // for longitudinal + grid measurement occasions

	array[Nevents] int cumu_events; // for longitudinal + grid measurement occasions (but for each event, rather than each individual)
	array[Nevents] int repme_events; // for longitudinal + grid measurement occasions  (but for each event, rather than each individual)

	matrix[Nlong, K] Y;

	vector[Nall] deltat; // gaps between all times for etas
  real delta_tx; // window defining how long tx effect lasts
  array[Nall] int repme_tx_times; // for each eta(t), how many prior tx were sent
  array[Nall] int cumu_tx_times; // for each eta(t), how many prior tx were sent
  array[Nall_tx] real time_a; // at what time were prior tx sent for each eta(t)
  array[Nall_tx] real lower_bound; // for each eta(t) and each prior tx, what is lower bound of integral in drift term
  array[Nall_tx] real upper_bound; // for each eta(t) and each prior tx, what is lower bound of integral in drift term
  array[Nall_tx] real time_t; // time of each eta
  
  vector[Nall] tau_factor; // number multiplied by tau to get tx effect (for haz and additive etas)
	
	vector[Nall] g_deltat_prev_event; // some function of time since most recent event

	array[Nevents] int status; // recurrent event status indicator
	
}

parameters {

  // measurement part
	vector<lower=0.000001>[K] lambda; // vector[K] lambda;
	real<lower=0.000001> sigma_lambda;
	vector<lower=0.000001>[K] sigma2_u;
	vector<lower=0.000001>[K] sigma2_e;

  // structural part
	matrix[Nall, P] eta;
	matrix[P, P] theta_ou;
	real<lower=-0.999999, upper=0.999999> rho; // correlation coefficient
	vector[P] tau; // intercept for tx-related drift

	// survival part
	real beta_0;
	real beta_1;
	real beta_2;
	real beta_3;
	real tau_tilde;
	
}

transformed parameters {

  // contstraint on theta_ou's eigenvalues
	real<lower=0.000001> constraint1;
	real<lower=0.000001> constraint2;

	constraint1 = theta_ou[1, 1] + theta_ou[2, 2];
	constraint2 = theta_ou[1, 1] * theta_ou[2, 2] - theta_ou[1, 2] * theta_ou[2, 1];

}

model{
  
  vector[Nall] haz; // Calculate hazard function at ALL times (grid pts + events)
  vector[Nevents] cumul_haz; // Calculate cumulative hazard function at event times
  vector[Nevents] haz_etimes; // Extract value of hazard function at event times
	
	matrix[P, P] V;
	
	// Prior
	lambda ~ normal(1, sigma_lambda);
	sigma_lambda ~ cauchy(0, 5);
	sigma2_u ~ cauchy(0, 5); // uses truncated half cauchu
	sigma2_e ~ cauchy(0, 5);

	to_vector(theta_ou) ~ normal(0, 10); //normal(0, 1)
	rho ~ uniform(-0.999999, 0.999999);
	to_vector(tau) ~ normal(0, 5); // tx effect on latent process

	beta_0 ~ normal(0, 5);
	beta_1 ~ normal(0, 5);
	beta_2 ~ normal(0, 5);
	beta_3 ~ normal(0, 5);
	tau_tilde ~ normal(0, 5); // tx effect on hazard (scalar)


  // Likelihood
	// At time=1
	V = [[1, rho], [rho, 1]];
	for (i in 1 : N){
		int k;
		k = cumu[i] - repme[i] + 1;
		eta[k] ~ multi_normal([0, 0]', V);
	}

	// Now is time = 2 to end
	for (i in 1 : N){
	  if (repme[i] > 1){
	    for (j in 2 : repme[i]){
			  int k;
			  int m;
			  vector[P] cond_mean;
			  // vector[P] drift_r;
			  matrix[P, P] drift_r = rep_matrix(0, P, P);
			  matrix[P, P] cond_covar;
			  matrix[P, P] term1_lower;
			  matrix[P, P] term2_lower;
			  matrix[P, P] term1_upper;
			  matrix[P, P] term2_upper;
			  k = cumu[i] - repme[i] + j;
			  for (l in 1 : repme_tx_times[k]){
			    m = cumu_tx_times[k] - repme_tx_times[k] + l;
		      // calculate intergral at lower bound
		      term1_lower = (1 - (lower_bound[m] - time_a[m]) / delta_tx) * matrix_exp(-1 * theta_ou * (time_t[m] - lower_bound[m])) / theta_ou;
          term2_lower = 1/delta_tx * matrix_exp(-theta_ou * (time_t[m] - lower_bound[m]));
          // calculate integral at upper bound
          term1_upper = (1 - (upper_bound[m] - time_a[m]) / delta_tx) * matrix_exp(-1 * theta_ou * (time_t[m] - upper_bound[m])) / theta_ou;
          term2_upper = 1/delta_tx * matrix_exp(-theta_ou * (time_t[m] - upper_bound[m]));
          // calculate addition to drift term for this specific tx time
          drift_r = drift_r + ((term1_upper + term2_upper) - (term1_lower + term2_lower));
			  }
			  cond_mean = matrix_exp(-deltat[k] * theta_ou) * to_vector(eta[k-1]);
			  cond_covar = V - matrix_exp(-deltat[k] * theta_ou) * V * matrix_exp(-deltat[k] * theta_ou');
			  eta[k] ~ multi_normal(cond_mean + drift_r * tau, cond_covar);
		  }
	  }
	}

	for (i in 1 : N){

	  int k;
	  int k_meas_start;
	  int k_meas_stop;
	  int ni = repme_meas[i];

	  vector[ni] eta1;
	  vector[ni] eta2;

	  matrix[ni, ni] Y1_var;
	  matrix[ni, ni] Y2_var;
	  matrix[ni, ni] Y3_var;
	  matrix[ni, ni] Y4_var;

	  k_meas_start = cumu_meas[i] - repme_meas[i] + 1;
	  k_meas_stop = cumu_meas[i];

	  eta1 = r_style_subset(eta, 1, meas_occ_rows[k_meas_start:k_meas_stop]);
	  eta2 = r_style_subset(eta, 2, meas_occ_rows[k_meas_start:k_meas_stop]);

	  Y1_var = add_diag(rep_matrix(sigma2_u[1], ni, ni), sigma2_e[1]);
	  Y2_var = add_diag(rep_matrix(sigma2_u[2], ni, ni), sigma2_e[2]);
	  Y3_var = add_diag(rep_matrix(sigma2_u[3], ni, ni), sigma2_e[3]);
	  Y4_var = add_diag(rep_matrix(sigma2_u[4], ni, ni), sigma2_e[4]);

    Y[k_meas_start:k_meas_stop, 1] ~ multi_normal(lambda[1]*eta1, Y1_var);
    Y[k_meas_start:k_meas_stop, 2] ~ multi_normal(lambda[2]*eta1, Y2_var);
    Y[k_meas_start:k_meas_stop, 3] ~ multi_normal(lambda[3]*eta2, Y3_var);
    Y[k_meas_start:k_meas_stop, 4] ~ multi_normal(lambda[4]*eta2, Y4_var);
	}


  // survival contribution to the likelhood
  for (j in 1 : Nall){ // Calc haz across grid pts
    haz[j] = exp(beta_0 + beta_1 * eta[j, 1] + beta_2 * eta[j, 2] + beta_3 * g_deltat_prev_event[j] + tau_tilde * tau_factor[j]);
    
  }
  
  for (i in 1 : Nevents){
    int k;
    cumul_haz[i] = 0;
    if (repme_events[i] > 1){
      for (j in 2 : repme_events[i]){
        k = cumu_events[i] - repme_events[i] + j;
        cumul_haz[i] = cumul_haz[i] + (haz[k] + haz[k-1])/2 * deltat[k];
      }
      haz_etimes[i] = haz[k]; // hazard at the event time
    } else {
      k = cumu_events[i] - repme_events[i] + 1;
      cumul_haz[i] = haz[k] * deltat[k]; 
      haz_etimes[i] = haz[k]; // hazard at the event time
    }
  }

  for (i in 1 : Nevents){
    target += status[i] * log(haz_etimes[i]) - cumul_haz[i];
  }

}

