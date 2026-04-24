################################################################################
# JOINT LONGITUDINAL RECURRENT EVENT MODEL WITH TX EFFECT
################################################################################

# this file references the programs needed to fit the model, calculate DIC and
# WAIC, and create posterior goodness-of-fit plots for the survival submodel
# using simulated data

# set working directory
my_wd <- ' '

# load useful functions (will be needed later)
source(paste0(my_wd, '/functions.R'))

################################################################################
# Which dataset to use?
################################################################################

# Number of individuals
N <- 100

# True parameter values (pick one setting)
setting <- 1
setting <- 2

# True treatment effect on the structural submodel (pick additive or drift)
tx_effect <- 'add' # eta(t) = oup + tx effect (additive impact)
tx_effect <- 'drift' # eta(t) has time-dependent drift (impact on derivative scale)

# True duration of treatment effect
delta_tx <- 0.5

# True hazard
haz_form <- 1 # h(t) = exp(beta0 + beta1 eta1 + beta2 eta2 + tilde_mu(t))
haz_form <- 2 # h(t) = exp(beta0 + beta1 eta1 + beta2 eta2 + beta3 time_since_prev_event + tilde_mu(t))

# True treatment effect as a function of time
mu_form <- 1 # exponential decay, tau * exp(-4 * time_since_tx)

################################################################################
# Which model to fit?
################################################################################

# (here we assume the fitted model is correctly specified, but it doesn't have to be)
mod_tx_effect <- tx_effect 
mod_haz_form <- haz_form

################################################################################
# Load the simulated data, fit the model, and save the results
################################################################################

# set seed (defines the simulated dataset that will be loaded, as well as stan)
g <- 1

# when fitting model model, need to specify the width of the grid used in the 
# cumulative hazard approximation
quad_gap <- 0.5

# load and set up data, then fit model
source(paste0(my_wd, '/fit_model/fit_ouf_rec.R'))

################################################################################
# Calculate DIC and WAIC
################################################################################

source(paste0(my_wd, '/check_fit/calc_dic_waic.R'))

################################################################################
# Generate goodness of fit plots
################################################################################

# need to define interval over which to generate posterior predictions of events
tstart <- 4
max_time <- max(long_dat$time)
tstop <- max_time # here, use max follow up as end of interval

source(paste0(my_wd, '/check_fit/ppc_survival.R'))






