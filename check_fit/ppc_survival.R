################################################################################
## POSTERIOR PREDICTIVE CHECK: EVENT PROBABILITIES VS OBSERVED EVENTS
################################################################################

library(mvtnorm)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(survival)
library(Rmpfr)
library(reshape2)
library(reReg)


source(paste0(my_wd, "functions.R"))


################################################################################
## Load posterior samples
################################################################################

fit <- readRDS(paste0(my_wd, 'fit_model/ouf_rec_s', setting,
                      '_TRUE_tx_effect_', tx_effect, '_haz_form_', haz_form,
                      '_FITTED_tx_effect_', mod_tx_effect, '_haz_form_', mod_haz_form,
                      '_g', g, '.RDS'))

# To calculate LPPD for individual i, we need to first load the fitted model
#  and set up the posterior samples of the parameters
P <- 2; K <- 4
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


# extract posterior parameter samples (assuming one chain)
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
## POSTERIOR EVENT PREDICTION
################################################################################

# number of posterior samples (subsample for speed)
max_s <- max(samples_all$iteration)
S_vec <- seq(1, max_s, by = 10)
S <- length(S_vec)

# we'll also need eta samples
post_etas <- as_draws_df(fit$draws('eta'))

# dataframes for saving posterior predictions: cumulative form and exact times
pp_results <- data.frame(id = c(), tstart = c(), tstop = c(),
                         pp_events = c(), s = c())
pp_results_times <- data.frame(id = c(), time = c(), s = c())


for (s in 1:S){
  cat('---------- s =', s, '----------\n')
  pp_output_s <- predict_events_i(data = dat,
                                  setting = setting,
                                  haz_form = mod_haz_form, s_id = s,
                                  tstart = tstart, tstop = tstop)
  # cumulative number of events within each grid window between tstart and tstop
  pp_events_s <- pp_output_s[[1]]
  pp_events_s$s <- s
  pp_results <- rbind(pp_results, pp_events_s)
  # events times of predictions between tstart and tstop (MORE IMPORTANT OUTPUT)
  pp_times_s <- pp_output_s[[2]]
  pp_times_s$s <- s
  pp_results_times <- rbind(pp_results_times, pp_times_s)
  # note that event times of -888 indicate ids not at risk between tstart and tstop
}


# For summaries, focus on actual event times and compare predicted events to
#  observed events between tstart and tstop using MCFs

# first, set up the observed data between tstart and tstop
surv_dat2 <- surv_dat %>%
    # censor at tstop
    mutate(event_status = ifelse(event_time < tstop, event_status, 0)) %>%
    mutate(event_time = ifelse(event_time < tstop, event_time, tstop)) %>%
    # select just time after tstart
    filter(event_time > tstart) %>%
    # and truncate time at tstop
    filter(event_time <= tstop)

# generate MCF
reObj_obs <- plot(with(surv_dat2,
                       Recur(time = event_time, id = id, event = event_status)),
                  mcf = TRUE)$data 


# now move on to the predicted event times
# also need to create censoring times for all at tstop
cens_dat3 <- data.frame(s = rep(1:S, each = dat$N), id = rep(1:dat$N, times = S))
cens_dat3 <- left_join(cens_dat3,
                       surv_dat %>% filter(event_status == 0) %>%
                         select(id, event_status, event_time),
                       by = 'id')
cens_dat3 <- cens_dat3 %>%
  mutate(event_censored = 1, event_time = tstop) %>%
  select(s, id, event_time, event_status, event_censored)
# and the combine censoring times and predicted event times
cumul_exp_events <- pp_results_times %>%
  mutate(event_time = time, event_status = 1, event_censored = 0) %>%
  select(s, id, event_time, event_status, event_censored) %>%
  rbind(cens_dat3 %>% filter(id %in% pp_results_times$id)) %>%
  arrange(s, id)
# # any id with an event time of -888 corresponds to a person who was not at risk between tstart and tstop,
# # we want to drop all observations (including censoring times for these people)
drop_ids <- unique((cumul_exp_events %>% filter(event_time == -888))$id)
cumul_exp_events <- cumul_exp_events %>%
  filter(!(id %in% drop_ids))


# now, calculate MCF for each s = 1, ..., S and combine across all s
# reObj_dat_preQ <- data.frame(time = c(), MCF = c(), s = c())
reObj_dat <- data.frame(time = c(), MCF = c(), s = c())

for (cur_s in 1:S){
  print(cur_s)
  reObj_obs_s <- with(cumul_exp_events %>% filter(s == cur_s),
                      Recur(time = event_time, id = id, event = event_status))
  reObj_dat_s <- plot(reObj_obs_s, mcf = TRUE)$data %>%
    select(time, MCF) %>%
    mutate(s = cur_s)
  reObj_dat <- rbind(reObj_dat, reObj_dat_s)
}
# if we group by time and calculate median and quantiles of MCF, then the resuling
#  step functions are very wiggly because groups have ~1 observation in each. to
#  calculate more reasonable looking median and quantiles, calculate rolling median
#  based on bins of 1 width (i.e., median of all MCF values within +/- 0.5 of current time)
rolling_window <- 1/max_time

reObj_dat_summ <- reObj_dat %>%
  arrange(time) %>%
  mutate(MCF_median = unlist(slider::slide_index(MCF, time, median, .before = rolling_window, .after = rolling_window)),
         MCF_lower = unlist(slider::slide_index(MCF, time, ~quantile(.x, prob = 0.05)[[1]], .before = rolling_window, .after = rolling_window)),
         MCF_upper = unlist(slider::slide_index(MCF, time, ~quantile(.x, prob = 0.95)[[1]], .before = rolling_window, .after = rolling_window)))
  # mutate(time = round(time, 3)) %>%
  # group_by(time) %>%
  # mutate(MCF_median = median(MCF),
  #        MCF_lower = quantile(MCF, probs = c(0.05)),
  #        MCF_upper = quantile(MCF, probs = c(0.95)))
  

plt <- ggplot() +
  coord_cartesian(xlim = c(0, tstop)) +
  geom_rect(alpha = 0.4, color = NA, fill = 'grey50',
            aes(xmin = 0, xmax = tstart, ymin = 0,
                ymax = max(reObj_dat_summ$MCF_upper, reObj_obs$MCF))) +
  geom_vline(aes(xintercept = (4*24)/max_time), linetype = 'dashed') +
  geom_ribbon(data = reObj_dat_summ,
              aes(x = time, ymin = MCF_lower, ymax = MCF_upper),
              fill = 'tomato2', alpha = 0.2) +
  geom_step(data = reObj_obs, aes(x = time, y = MCF), lwd = 1) +
  geom_step(data = reObj_dat_summ, aes(x = time, y = MCF_median),
            color = 'tomato2') +
  theme_bw() +
  labs(x = 'Time (days)', y = 'MCF')


# Goodness of fit plot
plt

# if you also want to save the data for the plot...
dat_to_save <- list(obs = reObj_obs,
                    pp = reObj_dat_summ)

saveRDS(dat_to_save,
        paste0(my_wd, 'ouf_rec_ppc_dat_v2_',
              'tx_effect_', mod_tx_effect, '_delta_tx_', delta_tx,
              '_haz', mod_haz_form, '_g', g, '_polysub_window', tstart, 'to', tstop, '.rds'))






