################################################################################
## CALCULATE MARGINAL LOG POSTERIOR PREDICTIVE DENSITY FOR USE WITH DIC, WAIC
################################################################################

library(mvtnorm)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(survival)


# distributions that we need to calculate are:
# log_lik_conditional = f_c(data | eta, params)
# eta_prior = p(eta | params)
# eta_post = p(eta | data)
# (functions for these distributions are defined in functions.R)

source(paste0(my_wd, "functions.R"))

################################################################################
## Load posterior samples
################################################################################

# samples that we need:
# all posterior draws
# etas, from the unconditional posterior density


# To calculate LPPD, we need to first load the fitted model
#  and set up the posterior samples of the parameters
fit <- readRDS(paste0(my_wd, 'data/fit_model/ouf_rec_s', setting,
                      '_TRUE_tx_effect_', tx_effect, '_haz_form_', haz_form,
                      '_FITTED_tx_effect_', mod_tx_effect, '_haz_form_', mod_haz_form,
                      '_g', g, '.RDS'))

# create vector of parameter names
P <- 2
if (mod_haz_form == 1){
  jm_params <- c('theta_ou[1,1]', 'theta_ou[2,1]',
                 'theta_ou[1,2]', 'theta_ou[2,2]', 'rho',
                 paste0('tau[', 1:P, ']'), 'tau_tilde',
                 'lambda[1]', 'lambda[2]', 'lambda[3]', 'lambda[4]',
                 'sigma2_u[1]', 'sigma2_u[2]', 'sigma2_u[3]', 'sigma2_u[4]',
                 'sigma2_e[1]', 'sigma2_e[2]', 'sigma2_e[3]', 'sigma2_e[4]',
                 'beta_0', 'beta_1', 'beta_2')
  
}else if (mod_haz_form == 2){
  jm_params <- c('theta_ou[1,1]', 'theta_ou[2,1]',
                 'theta_ou[1,2]', 'theta_ou[2,2]', 'rho',
                 paste0('tau[', 1:P, ']'), 'tau_tilde',
                 'lambda[1]', 'lambda[2]', 'lambda[3]', 'lambda[4]',
                 'sigma2_u[1]', 'sigma2_u[2]', 'sigma2_u[3]', 'sigma2_u[4]',
                 'sigma2_e[1]', 'sigma2_e[2]', 'sigma2_e[3]', 'sigma2_e[4]',
                 'beta_0', 'beta_1', 'beta_2', 'beta_3')
  
}


# extract posterior parameter samples (assuming one chain here)
raw_samples <- fit$draws(variables = c(jm_params))
samples <- matrix(raw_samples[,1,],
                  nrow = dim(raw_samples)[1],
                  ncol = dim(raw_samples)[3])
colnames(samples) <- jm_params
samples <- data.frame(samples)
samples$iteration <- 1:dim(raw_samples)[1]
samples_all <- melt(samples, c('iteration'))
colnames(samples_all) <- c('iteration', 'params', 'value')

################################################################################
## Calculate log f_m(y_i | theta^s)
################################################################################

# we need: log f_m(y_i | theta^s) and log f_m(y_i | mean(theta^s))

max_s <- 1000 # number of posterior samples
S_vec <- seq(1, max_s, by = 10) # subsample every 10th iteration for speed later
S <- length(S_vec)

# Then, for each individual i, we calculate the (approximate) unconditional
#. posterior distribution for eta_i
post_etas <- as_draws_df(fit$draws('eta'))
sampled_etas <- matrix(NA, nrow = 0, ncol = 5)
for (s in 1:max_s){
  post_etas2 <- post_etas[s,1:(dat$P*dat$Nall)]
  add_etas_mat <- data.frame(matrix(unlist(post_etas2), ncol = dat$P))
  post_samps <- cbind(id = all_times_template$id, s = s, 
                      time = all_times_template$time, add_etas_mat)
  colnames(post_samps) <- c('id', 's', 'time', 'eta1', 'eta2')
  sampled_etas <- rbind(sampled_etas, post_samps)
}
# calculate mean & variance
eta_i_means <- list(); eta_i_vars <- list()
for (i in 1:dat$N){
  print(i)
  eta_i <- sampled_etas %>% filter(id == i)
  temp_mean <- as.matrix(eta_i %>% group_by(time) %>%
                           summarise(eta1 = mean(eta1), eta2 = mean(eta2)) %>%
                           select(-c(time)))
  eta_i_wide <- reshape(eta_i %>% select(-c(id)), idvar = "s", timevar = "time",
                        direction = "wide") %>% select(-c(s))
  temp_covar <- cov(eta_i_wide)
  eta_i_means[[i]] <- as.vector(t(temp_mean))
  eta_i_vars[[i]] <- temp_covar
}

M <- 25 # for MC integration over latent process (bigger = better but SLOWER)

# ppdi <- matrix(NA, nrow = S, ncol = dat$N)
# function for calculating ppd_i
calc_ppdi <- function(s_id){
  s <- S_vec[s_id]
  cat('---- s =', s, '---- \n')
  # Set parameters to values of posterior sample s
  params_s <- list(
    'lambda' = samples_all[which(samples_all$iteration == s &
                                   grepl("lambda", samples_all$params)),]$value,
    'sigma2_u' = samples_all[which(samples_all$iteration == s &
                                     grepl("sigma2_u", samples_all$params)),]$value,
    'sigma2_e' = samples_all[which(samples_all$iteration == s &
                                     grepl("sigma2_e", samples_all$params)),]$value,
    'theta_ou' = matrix(samples_all[which(samples_all$iteration == s &
                                            grepl("theta_ou", samples_all$params)),]$value, 2, 2),
    'rho' = samples_all[which(samples_all$iteration == s &
                                grepl("rho", samples_all$params)),]$value,
    'tau' = samples_all[which(samples_all$iteration == s &
                                samples_all$params %in% c('tau.1.', 'tau.2.')),]$value,
    'beta_0' = samples_all[which(samples_all$iteration == s &
                                   samples_all$params == 'beta_0'),]$value,
    'beta_1' = samples_all[which(samples_all$iteration == s &
                                   samples_all$params == 'beta_1'),]$value,
    'beta_2' = samples_all[which(samples_all$iteration == s &
                                   samples_all$params == 'beta_2'),]$value,
    'tau_tilde' = samples_all[which(samples_all$iteration == s &
                                      samples_all$params == 'tau_tilde'),]$value
  )
  if (mod_haz_form == 2){
    params_s$beta_3 = samples_all[which(samples_all$iteration == s &
                                          samples_all$params == 'beta_3'),]$value
  }
  
  ## STEP 1: for a single set of posterior samples s, draw a list of length M that
  # contains the sampled eta values for each individual
  uncond_eta_M <- list()
  uncond_eta_m <- matrix(NA, nrow = 0, ncol = 2)
  for (m in 1:M){
    uncond_eta_m <- matrix(NA, nrow = 0, ncol = 2)
    for (i in 1:N){
      uncond_eta_m <- rbind(uncond_eta_m,
                            matrix(mvtnorm::rmvnorm(1, mean = eta_i_means[[i]],
                                           sigma = eta_i_vars[[i]]),
                                   ncol = 2, byrow = T))
    }
    colnames(uncond_eta_m) <- c('eta1', 'eta2')
    uncond_eta_M[[m]] <- uncond_eta_m
  }
  
  # STEP 2: Now that we have a set of M unconditional eta samples for each individual,
  # we need to calculate the values of the log densities for each m = 1, ..., M
  log_p_to_avg_over_m <- rep(0, dat$N)
  # optional: can set up this loop to run in parallel
  for (m in 1:M){
    cat('     m =', m, '\n')
    eta_m <- uncond_eta_M[[m]]
    
    log_lik_cond_m <- log_lik_conditional(data = dat, eta = eta_m,
                                          params = params_s, haz_form = mod_haz_form,
                                          surv_dat = surv_dat)
    
    log_eta_prior_m <- eta_prior(data = dat, eta = eta_m,
                                 params = params_s, tx_effect = mod_tx_effect)
    
    log_eta_post_m <- eta_post(data = dat, eta = eta_m,
                               eta_means = eta_i_means, eta_vars = eta_i_vars)
    
    # collect these log density values for this value of m:
    # lik_cond_m x p(eta_prior_m) / p(eta_post_m) = exp(log_lik_cond_m + log_eta_prior_m - log_eta_post_m)
    # but, let's keep this on the log scale for now...
    log_p_m <- log_lik_cond_m + log_eta_prior_m - log_eta_post_m
    log_p_to_avg_over_m <- log_p_to_avg_over_m + exp(log_p_m)
  }
  
  # save log_inside_term, since well need to sum exp(log_inside_term) across all
  #  s = 1, ..., S for each individual
  return(log_p_to_avg_over_m/M)
}

# for each set of posterior samples...(do this in parallel)
n_cores <- parallel::detectCores()
cat('using', n_cores, 'cores \n')
cl <- parallel::makeCluster(n_cores, outfile="")
doParallel::registerDoParallel(cl)

clusterExport(cl, c('samples_all', 'eta_i_means', 'eta_i_vars', 'data',
                    'surv_dat', 'mod_tx_effect', 'mod_haz_form'))
ppdi_list <- foreach(s_id = 1:S, .packages = c('dplyr', 'mvtnorm', 'Matrix')) %dopar% {
  calc_ppdi(s_id)
}
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)

# ppdi_list is a list of length S with each set of N llpi samples as elements
# want to convert this list into a S x N matrix
ppdi <- do.call(rbind,ppdi_list) 


################################################################################
## Calculate log f_m(y_i | mean(theta^s))
################################################################################

# calculate posterior mean for all parameters (not etas)
samples_all_mean <- samples_all %>%
  group_by(params) %>%
  summarise(value = mean(value))

# Set parameters to posterior mean
params_hat <- list(
  'lambda' = samples_all_mean[which(grepl("lambda", samples_all_mean$params)),]$value,
  'sigma2_u' = samples_all_mean[which(grepl("sigma2_u", samples_all_mean$params)),]$value,
  'sigma2_e' = samples_all_mean[which(grepl("sigma2_e", samples_all_mean$params)),]$value,
  'theta_ou' = matrix(samples_all_mean[which(grepl("theta_ou", samples_all_mean$params)),]$value, 2, 2),
  'rho' = samples_all_mean[which(grepl("rho", samples_all_mean$params)),]$value,
  'tau' = samples_all_mean[which(samples_all_mean$params %in% c('tau.1.', 'tau.2.')),]$value,
  'beta_0' = samples_all_mean[which(samples_all_mean$params == 'beta_0'),]$value,
  'beta_1' = samples_all_mean[which(samples_all_mean$params == 'beta_1'),]$value,
  'beta_2' = samples_all_mean[which(samples_all_mean$params == 'beta_2'),]$value,
  'tau_tilde' = samples_all_mean[which(samples_all_mean$params == 'tau_tilde'),]$value
)
if (mod_haz_form == 2){
  params_hat$beta_3 = samples_all_mean[which(samples_all_mean$params == 'beta_3'),]$value
}


## STEP 1: draw a list of length M that contains the sampled eta values for each individual
uncond_eta_M <- list()
uncond_eta_m <- matrix(NA, nrow = 0, ncol = 2)
for (m in 1:M){
  uncond_eta_m <- matrix(NA, nrow = 0, ncol = 2)
  for (i in 1:N){
    uncond_eta_m <- rbind(uncond_eta_m,
                          matrix(rmvnorm(1, mean = eta_i_means[[i]],
                                         sigma = eta_i_vars[[i]]),
                                 ncol = 2, byrow = T))
  }
  colnames(uncond_eta_m) <- c('eta1', 'eta2')
  uncond_eta_M[[m]] <- uncond_eta_m
}

# Now that we have a set of M unconditional eta samples for each individual,
# we need to calculate the values of the log densities for each m = 1, ..., M
# at theta = posterior mean (rather than sample s)
ppd_pt <- rep(0, dat$N)

# for each sample M, this can be done in parallel
n_cores <- parallel::detectCores()
cat('using', n_cores, 'cores \n')
cl <- parallel::makeCluster(n_cores, outfile="")
doParallel::registerDoParallel(cl)

clusterExport(cl, c('params_hat', 'eta_i_means', 'eta_i_vars', 'data',
                    'surv_dat', 'mod_tx_effect', 'mod_haz_form'))
ppd_list <- foreach(m = 1:M, .packages = c('dplyr', 'mvtnorm', 'Matrix')) %dopar% {
  cat('     m =', m, '\n')
  eta_m <- uncond_eta_M[[m]]
  log_lik_cond_m <- log_lik_conditional(data = dat, eta = eta_m,
                                        params = params_hat, haz_form = mod_haz_form,
                                        surv_dat = surv_dat)
  
  log_eta_prior_m <- eta_prior(data = dat, eta = eta_m,
                               params = params_hat, tx_effect = mod_tx_effect)
  
  log_eta_post_m <- eta_post(data = dat,
                             eta = eta_m,
                             eta_means = eta_i_means, eta_vars = eta_i_vars)
  
  # collect these log density values for this value of m:
  # lik_cond_m x p(eta_prior_m) / p(eta_post_m) = exp(log_lik_cond_m + log_eta_prior_m - log_eta_post_m)
  # but, let's keep this on the log scale for now...
  log_inside_term_m <- log_lik_cond_m + log_eta_prior_m - log_eta_post_m
  # log_inside_term_m <- mpfrArray(log_inside_term_m, prec = 64)
  # ppd_pt <- ppd_pt + exp(log_inside_term_m)/M
  exp(log_inside_term_m)/M
}
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)

# convert list to length-N vector
ppd_pt <- colSums(do.call(rbind, ppd_list))

################################################################################
## COMPILE RESULTS INTO LIST
################################################################################

res <- list('haz_form' = haz_form,
            'tx_effect' = tx_effect,
            'fitted_haz_form' = mod_haz_form,
            'fitted_tx_effect' = mod_tx_effect,
            'g' = g, M = M, S = S,
            lppd_i = log(ppdi), # dim = S x N (quantity 2)
            lppd_pt = log(ppd_pt) # dim = 1 x N (quantity 1)
            )

################################################################################
## CALCULATE DIC AND WAIC
################################################################################

# DIC = -2 * log f_m(y | theta_hat) + 2 p_DIC
llk_theta_hat <- sum(res$lppd_pt)
p_dic <- 2 * (sum(res$lppd_pt) - sum(res$lppd_i) / nrow(res$lppd_i))
dic <- -2 * llk_theta_hat + 2 * p_dic

# WAIC = - 2 lppd + 2 p_WAIC
lppd <- sum(log(colMeans(exp(res$lppd_i))))
# sample var
p_waic <- 0
for (i in 1:dat$N){
  a_si <- res$lppd_i[,i]
  a_s_bar <- mean(a_si)
  p_waic <- p_waic + sum((a_si - a_s_bar)^2) / (length(a_si) - 1)
}
waic <- -2 * lppd + 2 * p_waic

# SAVE RESULTS
res2 <- data.frame('haz_form' = haz_form,
                   'tx_effect' = tx_effect,
                   'fitted_haz_form' = mod_haz_form,
                   'fitted_tx_effect' = mod_tx_effect,
                   'g' = g, M = M, S = S,
                   waic = waic, dic = dic
)

write.csv(res2, row.names = FALSE,
          paste0(res_wd, 'lppd_results/info_crit_s', setting,
                 '_TRUE_tx_effect_', tx_effect, '_haz_form_', haz_form,
                 '_FITTED_tx_effect_', mod_tx_effect, '_haz_form_', mod_haz_form,
                 '_g', g, '.csv'))



