## predict.coxph.gammafrail (development version)

###  predict.coxph.gammafrail 1.0.0
- Version 1.0.0 provides the correct "baseline" (cumulative) hazard estimation
  since the `survival::survfit.coxph` function using the object `survival::coxph`
  with a  `survival::frailty()`term did not include in risk set the frailty term estimated
  by the Breslow estimator (the only one available for this function).
  Conditional predicted cumulative hazard and survival are consistent with `frailtyEM::predict.emfrail` using
   `frailtyEM::emfrail` for the shared gamma frailty semiparametric Cox models.




$$
\hat{H}_0(t) = \sum_{t_i \le t} \frac{d_i}{\sum_{j \in R_i} \exp\big(x_j^\top \hat{\beta} + \hat{\mu}_{c(j)}\big)}
$$

 $t_1 < t_2 < \dots < t_K\$ : unique event times  
- $d_i\$ : number of events at time $(t_i)$  
- $R_i =(j : T_j \ge t_i )$ : risk set at time $(t_i)$  
- $x_j$ : covariate vector for individual `j`  
- $\hat{\beta}$ : estimated fixed-effect coefficients  
- $\hat{u}_{c(j)}$ : estimated frailty for cluster $c(j)$  


  **NOTES:**
  - the current function version has not been (yet) neither implemented nor tested for left-truncation, non-linear terms (e.g., splines), stratification terms.
  - marginal estimated risk predictions are still experimental. We do not advice to use it yet.
  

###  predict.coxph.gammafrail 0.0.9000
- Adding NEWS file. The function provides the conditional and marginal predictions
  from a shared gamma frailty semiparametric Cox model using the function`survival::coxph`
  with a `survival::frailty()`term.
