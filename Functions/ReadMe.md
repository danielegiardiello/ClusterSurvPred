# Functions available 

+ rate_cens_bisection.R applies the iterative bisection procedure to define censor event times with a specified target proportion of censoring;

+ predict.coxph.gammafrail.R estimates conditional and marginal predictions at specific time horizons from the shared frailty gamma Cox proportional hazards models using
  the penalized partial likelihood through the function `survival::coxph`.
  The following [NEWS](https://github.com/danielegiardiello/ClusterSurvPred/blob/main/Functions/NEWS.md) file traces the updates of the function and provides important notes and disclaimers.

**CURRENT NOTES AND DISCLAIMERS**
1. The current predict.coxph.gammafrail.R function version has not been (yet) neither implemented nor tested for left-truncation, non-linear terms (e.g., splines), stratification terms.  
2. The marginal risk predictions estimation of predict.coxph.gammafrail.R is still experimental. We do not advice to use it yet.  

**Suggestions to improve these functions are kindly welcome**

