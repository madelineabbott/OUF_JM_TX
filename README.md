# Estimation of time-varying treatment effects in a joint model for longitudinal and recurrent event outcomes in mobile health data

This repository contains example code that carries out a model-based approach for estimating the effect of repeatedly delivered treatments in a micro-randomized trial (MRT) via an extension of a joint longitudinal-survival model.  The model and approach are described in ARXIVLINK.  This joint model, which summarizes multiple longitudinal outcomes as a smaller number of time-varying latent factors that capture the instantaneous risk of recurrent event outcomes, consists of three submodels: (i) a measurement submodel, (ii) a structural submodel, and (iii) a survival submodel.  We discuss different ways that these repeated treatment effects can be incorporated into the joint model; these different model specifications correspond to different mechanisms by which treatment is assumed to impact the longitudinal and event processes.  We take a Bayesian approach to inference in which we model the association between repeated treatments, multiple longitudinally measured outcomes, and recurrent events; inference is carried out via Stan (Carpenter et al., J. Stat. Soft., 2017).  We also present an approach for calculating information criteria that can be used for model selection.  The portion of code corresponding to longitudinal submodels (i) and (ii) builds on code written by Trung Dung Tran; the original code is available on the author's Github: https://github.com/tdt01/LOUmodels.

Two approaches to modeing the effect of repeated treatments:
+ additive (mean shift)
+ derivative scale (time-varying drift)

Approach to calculating information criteria:
+ DIC
+ WAIC

This repository contains simulated data for illustrative purposes in the [sim_dat](https://github.com/madelineabbott/OUF_JM_TX/tree/main/sim_dat) directory.

To fit the model, use either XXX (additive tx effect) or YYY (derivative tx effect).

To calculate information criteria, use XXX.

For questions, please contact mrabbott@umich.edu.
