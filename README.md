# Estimation of time-varying treatment effects in a joint model for longitudinal and recurrent event outcomes in mobile health data

This repository contains example code that carries out a model-based approach for estimating the effect of repeatedly delivered treatments in a micro-randomized trial (MRT) via an extension of a joint longitudinal-survival model.  The model and approach are described in ARXIVLINK.  This joint model, which summarizes multiple longitudinal outcomes as a smaller number of time-varying latent factors that capture the instantaneous risk of recurrent event outcomes, consists of three submodels: (i) a measurement submodel, (ii) a structural submodel, and (iii) a survival submodel.  We discuss different ways that these repeated treatment effects can be incorporated into the joint model; these different model specifications correspond to different mechanisms by which treatment is assumed to impact the longitudinal and event processes.  We take a Bayesian approach to inference in which we model the association between repeated treatments, multiple longitudinally measured outcomes, and recurrent events; inference is carried out via Stan (Carpenter et al., J. Stat. Soft., 2017). We also demonstrate how to calculate information criteria for model selection and present goodness-of-fit plots for assessing survival submodel calibration. The portion of code corresponding to longitudinal submodels (i) and (ii) builds on code written by Trung Dung Tran; the original code is available on the author's Github: https://github.com/tdt01/LOUmodels.

The repeated treatment effects can be incorporated into the structural submodel in two different ways:

+ We can assume that the impact of treatment is additive and directly shifts the latent process away from the constant mean for a short window of time after the treatment (called the **additive** approach).

+ We can model treatment as impacting the dynamics of the latent process through a time-varying drift term on the derivative scale (called the **drift** approach).



An overview of this repository is here:

1. Simulated data for illustrative purposes are available in the [sim_data](https://github.com/madelineabbott/OUF_JM_TX/tree/main/sim_data) directory.
2. To fit the model, use [fit_model/fit_ouf_rec.R](https://github.com/madelineabbott/OUF_JM_TX/blob/main/fit_model/fit_ouf_rec.R).
3. After fitting the model, models can be compared using information criteria, with code provided in [check_fit/calc_dic_waic.R](https://github.com/madelineabbott/OUF_JM_TX/blob/main/sim_data/check_fit/calc_dic_waic.R).
4. Goodness-of-fit plots for assessing the survival submodel calibration can be generated using code in [check_fit/ppc_survival.R](https://github.com/madelineabbott/OUF_JM_TX/blob/main/sim_data/check_fit/ppc_survival.R).


For questions, please contact mabbott@sdac.harvard.edu.
