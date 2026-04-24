################################################################################
# FIT JOINT MODEL WITH LATENT PROCESS AND TREATMENT EFFECT
################################################################################

library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(survival)

# set seed
set.seed(22210 + g)

# # NUMBER OF INDIVIDUALS
# N <- 100
# 
# # TRUE PARAMETER VALUES
# setting <- 1
# setting <- 2
# 
# # GRID WIDTH USED TO APPROXIMATE CUMULATIVE HAZARD
# quad_gap <- 0.5
# # if quad_gap == 999, then only measurement occasions will be used (no additional grid points)
# 
# # WHERE DOES TX IMPACT THE STRUCTURAL SUBMODEL (ADDITIVE OR DRIFT)
# tx_effect <- 'add' # eta(t) = oup + tx effect
# tx_effect <- 'drift' # eta(t) has time-dependent drift
# 
# ## WHAT DOES THE HAZARD FUNCTION LOOK LIKE?
# haz_form <- 1 # h(t) = exp(beta0 + beta1 eta1 + beta2 eta2 + tilde_mu(t))
# haz_form <- 2 # h(t) = exp(beta0 + beta1 eta1 + beta2 eta2 + beta3 time_since_prev_event + tilde_mu(t))
# 
# # HOW DOES TREATMENT EFFECT CHANGE OVER TIME
# mu_form <- 1 # exponential decay, tau * exp(-4 * time_since_tx)
# 
# # WHAT IS DURATION OF TREATMENT EFFECT?
# delta_tx <- 0.5
# 
# ### DEFINE MODEL TO BE FIT ###
# # (assuming fitted model is correctly specified here)
# mod_tx_effect <- tx_effect
# mod_haz_form <- haz_form

source(paste0(my_wd, "functions.R"))

################################################################################
# LOAD DATA
################################################################################

tx_dat <- read.csv(paste0(my_wd, 'sim_data/tx_dat_s', setting,
                          '_tx_effect_', tx_effect, '_N', N, '_g', g, '.csv'))

long_dat <- read.csv(paste0(my_wd, 'sim_data/long_dat_s', setting,
                            '_tx_effect_', tx_effect, '_N', N, '_g', g, '.csv'))

surv_dat <- read.csv(paste0(my_wd, 'sim_data/surv_dat_s', setting,
                           '_tx_effect_', tx_effect,'_haz_form_', haz_form,
                           '_N', N, '_g', g, '.csv'))

################################################################################
# DEFINE GRID FOR SURVIVAL FUNCTION APPROXIMATION
################################################################################

t_max <- max(surv_dat$event_time)
num_quad <- ceiling(t_max / quad_gap)
N <- length(unique(surv_dat$id))

# these will be the extra grid points used when calculating the cumulative hazard
aug_dat <- data.frame(id = rep(1:N, each = num_quad + 1),
                      time = rep(seq(0, t_max, length.out = num_quad + 1),
                                 times = N),
                      time_type = 'quad')

# these are the grid points defined by the measurement occasions
long_dat_temp <- long_dat %>%
  dplyr::select(id, time) %>%
  mutate(time_type = 'meas')

# and we include a grid point at the event times
surv_dat_temp <- surv_dat %>%
  mutate(time = event_time, time_type = 'event') %>%
  dplyr::select(id, time, time_type)


if (quad_gap == 999){
  # no grid, just measurement occasions and event time
  all_times_template <- rbind(long_dat_temp, surv_dat_temp) %>%
    arrange(id, time, time_type)
}else{
  # otherwise, combine timepoints defined by extra grid points, measurement
  # times, and event times for use approximating the cumulative hazard later
  all_times_template <- rbind(aug_dat, long_dat_temp, surv_dat_temp) %>%
    arrange(id, time, time_type)
}


# remove grid points after survival time
final_event_time <- surv_dat %>%
  group_by(id) %>%
  arrange(event_time, .by_group = T) %>%
  slice(n()) %>%
  ungroup()
all_times_template <- left_join(all_times_template, final_event_time, by = 'id')

all_times_template <- all_times_template %>%
  filter(time <= event_time) %>%
  dplyr::select(id, time, time_type)

# and filter out some other unnecessary grid points
all_times_template <- all_times_template %>%
  group_by(id) %>%
  mutate(keep1 = ifelse(time_type %in% c('meas', 'event'), 1, 0),
         keep2 = ifelse((time_type == 'quad' & lag(time_type) != 'quad' &
                           abs(time - lag(time)) < 1e-8), 0, 1),
         keep3 = ifelse((time_type == 'quad' & lead(time_type) != 'quad' &
                           abs(lead(time) - time) < 1e-8), 0, 1),
         keep4 = ifelse((time_type == 'event' & lead(time_type) == 'quad' &
                           abs(time - lead(time)) < 1e-8), 1, 0)) %>%
  mutate(keep1 = ifelse(is.na(keep1), 0, keep1),
         keep2 = ifelse(is.na(keep2), 0, keep2),
         keep3 = ifelse(is.na(keep3), 0, keep3),
         keep4 = ifelse(is.na(keep4), 0, keep4)) %>%
  mutate(keep = as.numeric((keep1 + keep2 * keep3 + keep4) > 0)) %>%
  ungroup() %>%
  filter(keep == 1)

# if a measurement and event happen at the same time, then create new time_type
all_times_template <- all_times_template %>%
  mutate(time_type_old = time_type) %>%
  arrange(id, time, time_type) %>%
  group_by(id) %>%
  mutate(time_type = case_when((time_type == 'meas' &
                                  lag(time_type) == 'event' &
                                  (time - lag(time)) < 1e-8) ~ 'meas_event',
                               (time_type == 'event' &
                                  lead(time_type) == 'meas' &
                                  (lead(time) - time) < 1e-8) ~ 'meas_event',
                               TRUE ~ time_type_old)) %>%
  ungroup() %>%
  distinct(id, time, time_type, .keep_all = TRUE)

# if quad. points are too close to meas. points, then drop the quad. point
# note that we define "too close" as within 30% of quad_gap
nrow_orig_template <- nrow(all_times_template)
all_times_template <- all_times_template %>%
  group_by(id) %>%
  mutate(min_gap = abs(pmin(time - lag(time, default = -99),
                            time - lead(time, default = -99)))) %>%
  ungroup() %>%
  mutate(drop = ifelse(time_type == 'quad' & min_gap <= 0.3*quad_gap, 1, 0)) %>%
  filter(drop == 0) %>%
  dplyr::select(-c(min_gap, drop))


# create variables that contains gaps between cleaned-up grid points
all_times_template <- all_times_template %>%
  group_by(id) %>%
  mutate(deltat = time - lag(time, default = 0)) %>%
  ungroup()

# add in event-specific IDs
all_times_template <- all_times_template %>%
  mutate(event_id = 1 + cumsum(time_type %in% c('event', 'meas_event'))) %>%
  group_by(event_id) %>%
  mutate(deltat_event = time - lag(time, default = 0)) %>%
  ungroup() %>%
  mutate(event_id = ifelse(time_type %in% c('event', 'meas_event'),
                           event_id - 1, event_id))



################################################################################
# DEFINE INTERVENTION TIMES & PRECALCULATE FACTOR CAPTURING IMPACT ON LATENT PROCESS
################################################################################

# if tx intervals are allowed to overlap, then we need to precalculate the
#  term onto which tau is multiplied to make mu

# IF TX IS MODELED AS HAVING AN ADDITIVE IMPACT ON THE LATENT PROCESS...
# if tx intervals are allowed to overlap, then we need to precalculate the
#  term onto which tau is multiplied to make mu.  That is, we precalculate:
# tau_factor = (1 - (time-tx_time)/drift_window) but account for lap in tx windows
tau_factor <- c() # this factor is the same for eta and hazard, but the taus are different
for (cur_id in 1:length(unique(all_times_template$id))){
  long_dat_i <- all_times_template %>% filter(id == cur_id)
  tx_dat_i <- tx_dat %>% filter(id == cur_id)
  for (cur_time in long_dat_i$time){
    tau_factor <- c(tau_factor,
                    calculate_tau_factor(cur_time = cur_time,
                                         tx_time = tx_dat_i$tx_time,
                                         delta_tx = delta_tx))
  }
}
all_times_template$tau_factor = tau_factor



# IF TX IS MODELED VIA DRIFT (i.e. ASSUMED TO IMPACT THE DERIVATIVE OF THE LATENT PROCESS)...
time_a_vec <- c(); lower_bound_vec <- c();
upper_bound_vec <- c(); time_t_vec <- c()
repeated_eta_id_vec <- c() # id specific to each eta_i(t)
repeated_eta_id <- 0
for (cur_id in 1:length(unique(all_times_template$id))){
  long_dat_i <- all_times_template %>% filter(id == cur_id)
  tx_dat_i <- tx_dat %>% filter(id == cur_id)
  # assume drift is 0 at t0 (j = 1)
  time_t_vec <- c(time_t_vec, 0)
  time_a_vec <- c(time_a_vec, 0)
  lower_bound_vec = c(lower_bound_vec, 0)
  upper_bound_vec = c(upper_bound_vec, 0)
  repeated_eta_id <- repeated_eta_id + 1
  repeated_eta_id_vec <- c(repeated_eta_id_vec, repeated_eta_id)
  # then start j index at 2
  for (j in 2:nrow(long_dat_i)){
    cur_time <- long_dat_i$time[j]
    cur_prior_tx <- tx_dat_i %>% filter(long_dat_i$time[j-1] - delta_tx <= tx_time  & tx_time < cur_time)
    # if no tx were sent recently, then we don't need to calculate drift since drift = 0
    if (nrow(cur_prior_tx) == 0){
      time_t_vec <- c(time_t_vec, 0)
      time_a_vec <- c(time_a_vec, 0)
      lower_bound_vec = c(lower_bound_vec, 0)
      upper_bound_vec = c(upper_bound_vec, 0)
      repeated_eta_id <- repeated_eta_id + 1
      repeated_eta_id_vec <- c(repeated_eta_id_vec, repeated_eta_id)
    }else if (nrow(cur_prior_tx) > 0){
      repeated_eta_id <- repeated_eta_id + 1
      for (k in 1:nrow(cur_prior_tx)){
        # for each recently delivered tx, need tx time and upper/lower integral bounds
        time_a <- cur_prior_tx$tx_time[k]
        time_s <- long_dat_i$time[j-1]
        time_t_vec <- c(time_t_vec, cur_time)
        time_a_vec <- c(time_a_vec, time_a)
        lower_bound_vec = c(lower_bound_vec, max(time_a, time_s))
        upper_bound_vec = c(upper_bound_vec, min(time_a + delta_tx, cur_time))
        # record id that is unique to each eta_i(t) (tells us how many recent txs we need to sum over in drift term)
        repeated_eta_id_vec <- c(repeated_eta_id_vec, repeated_eta_id)
      }
    }
  }
}

# total number of recent tx for each eta(t)
repme_tx_times <- table(repeated_eta_id_vec) 
# cumulative number of recent tx for each eta(t)
cumu_tx_times <- cumsum(table(repeated_eta_id_vec)) 



################################################################################
# PRECALCULATE FACTOR CAPTURING IMPACT OF INTERVENTION ON SURVIVAL SUBMODEL
################################################################################

# if hazard depends on some function of time since most recent event (e.g., logit),
# then we need to precalculate this value
g_deltat_prev_event <- c()
# hazard depends also on time since previous event
for (cur_id in 1:N){
  all_times_template_i <- all_times_template %>% filter(id == cur_id)
  surv_dat_i <- surv_dat %>% filter(id == cur_id)
  most_recent_event_time <- -999
  for (j in 1:nrow(all_times_template_i)){
    # determine time of most recent event
    cur_time <- all_times_template_i$time[j]
    prior_event_time <- surv_dat_i %>%
      filter(event_time < cur_time & event_status == 1) %>%
      slice(n())
    no_prior_events <- (nrow(prior_event_time) == 0)
    if (no_prior_events == FALSE){
      most_recent_event_time <- prior_event_time$event_time
    }
    time_since_prev_event <- cur_time - most_recent_event_time
    if (setting == 1){
      g_time_since_prev_event <- 1 / (1 + exp(4*(time_since_prev_event - 2)))
    }else if (setting == 2){
      g_time_since_prev_event <- 1 / (1 + exp(1.5*(time_since_prev_event - 2)))
    }
    g_deltat_prev_event <- c(g_deltat_prev_event, g_time_since_prev_event)
  }
}
all_times_template$g_deltat_prev_event = g_deltat_prev_event



################################################################################
# SET UP LIST OF DATA FOR STAN
################################################################################

dat = list(
  # number of meas. occasions + grid points
  Nall = nrow(all_times_template),
  # number of meas. occasions only (no grid points)
  Nlong = sum(all_times_template$time_type %in% c('meas', 'meas_event')),
  # number of subjects
  N = length(unique(long_dat$id)),
  # number of events
  Nevents = nrow(surv_dat),
  # number of items
  K = 4,
  # number of latent factors
  P = 2,
  # ID column from long. data (needed for random intercepts)
  measID = all_times_template$id[which(all_times_template$time_type %in% c('meas', 'meas_event'))],
  # cumulative # of meas occ. + grid points per person (for generating etas)
  cumu = cumsum(c(table(all_times_template$id))),
  # total # meas occ. + grid points per person (for generating etas)
  repme = c(table(all_times_template$id)), 
  
  # row ID corresponding to etas at measurement occasions of Y
  meas_occ_rows = c(1:nrow(all_times_template))[which(all_times_template$time_type %in% c('meas', 'meas_event'))],
  cumu_meas = cumsum(c(table(all_times_template$id[which(
    all_times_template$time_type %in% c('meas', 'meas_event'))]))),
  # total # meas occ. per person (for Ys)
  repme_meas = c(table(all_times_template$id[which(
    all_times_template$time_type %in% c('meas', 'meas_event'))])), 
  
  # for recurrent events
  cumu_events = cumsum(c(table(all_times_template$event_id))),
  repme_events = c(table(all_times_template$event_id)),
  
  # row #s corresponding to meas occ among meas occ + grid points (for obs. outcome)
  meas_occ = c(1:nrow(all_times_template))[which(all_times_template$time_type %in% c('meas', 'meas_event'))],

  # for treatment effect in structural submodel (if modeled via additive)
  tau_factor = all_times_template$tau_factor,
  # for treatment effect in structural submodel (if modeled via drift)
  time_t = time_t_vec, # time of each eta, repeated for tx
  delta_tx = delta_tx, # how long tx effect lasts
  Nall_tx = length(time_a_vec), # total num of tx used in drift term, including repeats
  repme_tx_times = repme_tx_times, 
  cumu_tx_times = cumu_tx_times,
  time_a = time_a_vec, # for each eta(t), when was each prior tx sent
  lower_bound = lower_bound_vec, # for each eta(t) and each prior tx, what is lower bound of integral in drift term
  upper_bound = upper_bound_vec, # for each eta(t) and each prior tx, what is upper bound of integral in drift term            
  
  # observed long. outcome
  Y = long_dat[,paste0('y', 1:4)],
  # event indicator
  status = surv_dat$event_status,
  # gap times (between full grid formed by meas. occ. + grid points)
  deltat = all_times_template$deltat,
  # hazard depends on some function of time since prev. event
  g_deltat_prev_event = all_times_template$g_deltat_prev_event,
  # times of grid points for calc of baseline haz (if time-varying)
  cur_time = all_times_template$time + 0.0001)

# initialize parameters (alternatively, could use two-stage approach)
init_params <- list(
  'lambda' = rep(1, 4),
  'sigma2_u' = rep(0.5, 4),
  'sigma2_e' = rep(0.5, 4),
  'theta_ou' = array(c(1, 0.5, 0.5, 1), dim = c(2, 2)),
  'rho' = -0.5,
  'tau' = c(1, -1),
  'beta_0' = -1, 'beta_1' = -1, 'beta_2' = 1, 'beta_3' = 1,
  'tau_tilde' = -0.5
)

################################################################################
# FIT MODEL
################################################################################

# load and compile model
file <- file.path(paste0(my_wd, 'fit_model/fit_OUF_REC_tx_effect_', mod_tx_effect,
                         '_haz_form_', mod_haz_form, '.stan'))
mod <- cmdstan_model(file)

# Run MCMC using the 'sample' method
mcmc <- mod$sample(
  data = dat,
  seed = g + 380, 
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 1000, # default = 1000
  iter_sampling = 1000, # default = 1000
  init = list(init_params),
  save_warmup = F # save warm up iterations
)

################################################################################
# SAVE SAMPLES
################################################################################

mcmc$save_object(file = paste0(my_wd, 'fit_model/ouf_rec_s', setting,
             '_TRUE_tx_effect_', tx_effect, '_haz_form_', haz_form,
             '_FITTED_tx_effect_', mod_tx_effect, '_haz_form_', mod_haz_form,
             '_g', g, '.RDS'))

