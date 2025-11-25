
#' @description
#' 
#' Function to estimate conditional and marginal 
#' survival probabilities at specified time horizons
#' for shared Gamma frailty Cox proportional hazards models
#' 
#' @author Daniele Giardiello
#' @param model cox model using survival::coxph including frailty gamma term
#' The event of interest must have the indicator variable equal to 1.
#' The model survival::coxph() must specify y = T survival::coxph(..., y = T)
#' The user must define a reference value
#' for all predictors, especially the continuous predictors
#' Example: age50 = age - 50 or age_c = age - mean(age)
#' Then "the baseline" corresponds to a subject with 
#' the specified reference predictors values. 
#' Ties/methods must be "breslow". 
#' @param newdata newdata (for validation). This should have the same variables of the model.
#' Predictors with levels should be coded with exactly the same way and levels of the predictors
#' with factors used to develop the model
#' @param cluster cluster variable. Only one variable possible
#' @param times prediction horizon time(s)
#' @examples
#' # example code
#' data(lung, package = "survival")
#' lung$age_c <- lung$age - mean(lung$age)
#' cox_frailty_lung <- survival::coxph(Surv(time, status) ~ age_c + frailty(inst,
#' distribution = "gamma"), 
#' data = lung,
#' x = TRUE,
#' y = TRUE,
#' model = TRUE,
#' ties = "breslow")
#' predict.coxph.gammafrail(model = cox_frailty_lung,
#' newdata = lung_sel,
#' cluster = "inst",
#' times = c(100, 200, 300))
#' @references Collett D - Modelling Survival data in Medical Research, 4th edition chapter 10 (formula 10.9 and 10.12)


# **Disclaimers**
# 1. The current functions provided the conditional estimated
# predictions using the shared gamma frailty Cox model using
# survival::coxph.
# Left-truncation, non-linear terms and stratification terms
# have not been implemented and tested yet. 

# 2. Marginal prediction estimation are still
# experimental. We do not advice to use it yet.
# Your suggestions to improve the function are kindly and warmly welcome. 


predict.coxph.gammafrail <- function(model,
                                     newdata,
                                     cluster,
                                     times) {
  
  # sanity check message warning
  on.exit(base::message("
  Sanity check warning: 
  please center all variables before fitting the model or assign reference values
                for continuous variable. 
                For example age50 = age - 50.
                newdata should contain the centered/referenced values"))
  
  # Stopping messages
  if (!inherits(model, "coxph"))
    stop("`model` must be a survival::coxph fit. 
         The model must specify y = T and Breslow method. 
         Please center all variables before fitting the model or assign reference values
         for continuous variable. For example age50 = age - 50")
  
  if (missing(cluster))
    stop("Please supply the name of the cluster variable (character).")
  
  if (!is.character(cluster) || length(cluster) != 1L)
    stop("`cluster` must be a single character string giving the cluster column name.")
  
  if (!cluster %in% names(newdata))
    stop("`cluster` not found in `newdata`.")
  
  
  # Start function ----
  # Step0: 
  # Formula with all terms including cluster as it is a fixed covariate
  frm_allterms <- reformulate(base::all.vars(model$formula[[3]]))
  
  # Select the complete cases of fixed + frailty terms
  # of the development data
  # (same data used to develop the model)
  # The argument survival::coxph(..., x = T) provides
  # the development data
  # NOTE: the development data is necessary to 
  # estimate the "baseline" (cumulative) hazard H0 later
  
  devdata <- stats::model.frame(frm_allterms,
                                data = base::eval(model$call$data, 
                                                  envir = environment(formula(model))), 
                                # NOTE: data.frame(model$x) did not
                                # provide the same variable names
                                # (e.g. frailty(id_cluster) instead of id_cluster)
                                # Alternative: add devdata as a new argument of the
                                # function. 
                                # The alternative option is probably the most
                                # intuitive 
                                na.action = na.omit) # not sure if necessary but safer
  
  
  # Select only the complete cases of fixed + frailty terms
  # of the newdata
  
  newdata <- stats::model.frame(frm_allterms,
                                data = newdata,
                                na.action = na.omit)
  
  # Step 1: survfit for
  # a baseline cumulative hazard / survival
  
  # Survfit of the coxph with frailty
  # sf_model <- survival::survfit(model)
  
  # NOTE: the survival::survfit with frailty does not accept newdata
  # thus, it is important to center all variables or give them a reference.
  # In this setting survfit$cumhaz will directly represent our baseline 

  # NOTE: [20251029]: survfit() with frailty with no newdata 
  # seems not to include the frailty terms in the "baseline"
  # cumulative hazard estimation for a reference subject
  
  # Step 1: linear predictors of fixed effects
  
  # Estimate linear predictors only
  # for fixed-effect effects without frailties
  # This will be also useful later when 
  # marginal predictions are estimated

  
  # Since with clusters size < vs >= 5
  # frailty terms are treated as coefficients or not,
  # I prefer calculate lp for the only fixed terms only manually
  
  # lp calculation for marginal
  # save all model coefficients
  # for cluster < 5, frailty terms are included in coefficients (not practical)
  coefs <- model$coefficients
  
  # Remove frailty named "gamma"
  # NOTE/WARNING: 
  # if a variables called gamma it must be a problem
  # to be adjusted
  coefs_fixed <- coefs[!base::grepl("gamma", names(coefs))]
  
  # create formula of fixed terms only
  terms_vec <- attr(terms(model$formula), 
                    "term.labels")
  
  # Remove specials using grep/grepl
  # NOTE: working in progress
  specials <- c("frailty", "strata", "cluster", "offset", "tt")
  fixed_terms <- terms_vec[!grepl(paste(specials, collapse = "|"), terms_vec)]

  
  formula_fixed <- stats::reformulate(fixed_terms,
                                      response = NULL)
  
  # create design matrix of fixed term only
  # of the new data
  X_fixed <- model.matrix(formula_fixed,
                          newdata)

  X_fixed <- X_fixed[, colnames(X_fixed) != "(Intercept)"]

  # estimate linear predictor X*beta fixed term only
  lp_fixed <- X_fixed %*% cbind(coefs_fixed)
  
  # this should be equivalent to this (when cluster > 5)
  # lp_nofrail <- cbind(predict(model,
  #                             newdata = newdata,
  #                             type = "lp"))
  
  # create design matrix of fixed term only
  # of the development data
  # (data used to develop the model)
  # useful to estimate the "baseline" cumulative hazard
  X_fixed_dev <- model.matrix(formula_fixed,
                              devdata)
  
  X_fixed_dev <- X_fixed_dev[, colnames(X_fixed_dev) != "(Intercept)"]
  
  # estimate linear predictor X*beta fixed term only
  lp_fixed_dev <- X_fixed_dev %*% cbind(coefs_fixed)
  
  
  # Step 2: linear predictors for frailty terms
  
  # NOTE/WARNING:
  # for cluster < 5 predict.coxph(model, newdata) includes also frailty
  # and frailty terms must not be included for marginal predictions.
  # This leads to wrong estimations of marginal predictions 
  # and the corresponding performance metrics in riskRegression::Score()
  
  # For frailty with not model$frail (where clusters < 5)
  # add the model$frail and keep only fixed-effect coefficients
  if(is.null(model$frail)) {
    coefs_model <- model$coefficients
    model$frail <- unname(coefs_model[base::grep("gamma", names(coefs_model))])
  }
  
  # Force cluster as a factor in the newdata
  # number of clusters in the newdata 
  
  # otherwise a factor
  # NOTE: for newdata with one cluster, factor does not work
  # newdata[[cluster]] <- base::as.character(newdata[[cluster]])
  
  # newdata 
   newdata[[cluster]] <- base::factor(newdata[[cluster]],
                                      levels = unique(newdata[[cluster]]))
   
  # development data
   devdata[[cluster]] <- base::factor(devdata[[cluster]],
                                      levels = unique(devdata[[cluster]]))
  
  # newdata[[cluster]] <- base::factor(newdata[[cluster]],
  #                                    levels = 1:length(model$frail))
   
  # Create the design matrix for frailty terms
  all_terms <- base::all.vars(model$formula[[3]])
  frailty_terms <- c(-1, all_terms[all_terms == cluster])
  
  # design matrix for the frailty terms
  # of the new data
  X_frail <- stats::model.matrix(stats::reformulate(frailty_terms),
                                 response = NULL,
                                 data = newdata)
  
  # design matrix for the frailty terms
  # of the development data (data used to develop the model)
  # useful to estimate the "baseline" cumulative hazard H0
  X_frail_dev <- stats::model.matrix(stats::reformulate(frailty_terms),
                                     response = NULL,
                                     data = devdata)
  
  # Frailty terms linear predictors
  # new data
  lp_frail <- X_frail %*% cbind(model$frail)
  
  # Frailty terms linear predictors
  # development data
  lp_frail_dev <- X_frail_dev %*% cbind(model$frail)
  
  # Total linear predictor
  # new data
  lp_total <- lp_fixed + lp_frail
  
  # Total linear predictor
  # development data
  lp_total_dev <- lp_fixed_dev + lp_frail_dev
  
  # Step 3: 
  # estimate the "baseline" cumulative hazard
  # for a reference subject 
  # NOTE: we use the development data
  # NOTE: recalibration might be possible using
  # the linear predictors (lp) of the development data
  # but the risk setting of the newdata (i.e., validation)
  # It might be considered as future developments
  
  risk_dev <- exp(lp_total_dev) # exp(lp + frailty)

  # Unique event times of the
  # model fitted using the development data
  
  model_times <- model$y[, "time"] # event time for development data
  model_status <- model$y[, "status"] # event indicator for development data
  
  # NOTE: the object model must specify y = T
  # Example: survival::coxph(..., y = T)
  # NOTE: status == 1 must be the event of reference
  
  # Extract the unique event time
  # from the development data
  # used to develop the model
  unique_event_times <- sort(unique(model_times[model_status == 1]))
  n_times <- length(unique_event_times)
  
  # Initialize vectors
  h0 <- numeric(n_times)    # incremental baseline hazard
  d_i <- numeric(n_times)   # number of events at each time
  
  for (i in seq_along(unique_event_times)) {
    t <- unique_event_times[i]
    d_i[i] <- sum(model_status[model_times == t])
    
    # Risk set: all individuals still at risk at time t
    risk_set <- which(model_times >= t)
    
    # Breslow baseline hazard
    h0[i] <- d_i[i] / sum(risk_dev[risk_set])
  }
 
  # Cumulative "baseline" hazard
  # Create the sf model
  H0 <- cumsum(h0)
  
  # Cumulative hazards
  # at fixed time horizons
  # H_0(t)
  cumhaz0_thor <- stats::approx(
    x = unique_event_times,
    y = H0,
    yleft = 0, # y values below min(x). Cumulative hazard zero if below min(time)
    yright = NA, # y values above max(x). Cumulative hazard NA if above max(time)
    xout = times,
    method = "constant",
    f = 0,
    rule = 1)$y
  
  # NOTE: 
  # stats:approx should work fine since
  # all values from survival::survfit are unique
  # for the event time of development data
  
  # H(t) using the newdata (i.e., validation data)
  cumhaz_thor <- exp(lp_total) %*% cumhaz0_thor
  
  # conditional S(t)
  S_cond_t <- exp(-cumhaz_thor)
  
  # marginal S(t)
  # NOTE: marginal predictions needs
  # to be further investigated
  # to check if they are corrected
  # 
  # NOTE/QUESTION
  # does the marginal use the baseline hazard including the frailty terms?
  
  theta_hat <- model$history[[1]]$theta
  S_marg_t <- (1 + theta_hat * (exp(lp_fixed) %*% cumhaz0_thor))^(-1/theta_hat)
  
  # Output 
  # NOTE: to be better defined
  res <- list("conditional" = S_cond_t,
              "marginal" = S_marg_t)
  
  return(res)
  # Report only conditional predictions
  # return(S_cond_t)
  
}



