## 28th July 2025
# Goal: simulating event times with a specified proportion of censoring 
# Authors: D Giardiello E Ratti

# Attempt to create a function
#' Function to simulate event times with a specified proportion of censoring
#' 
#' @param n number of observations 
#' @param tolerance tolerance for the iterative bisection procedure
#' @param prop.cens.target target censoring proportion
#' @param prop.cens.empirical empirical censoring proportion. It must be different than the target proportion
#' @param int_low
#' @param int_high
#' @param event_time simulated event times without censoring (for example using Cox-Weibull). 
#' The length of the vector must be equal to n argument (number of observations)
#' @param max_iter max number of iterations to avoid long iterations and computations. Default 10000
#' @param dgt number of digits to show the output. Default is three digits
#' @param seed set.seed

rate_cens_bisection <- function(
    n,
    tolerance = 0.001,
    prop.cens.target,
    prop.cens.empirical,
    int_low,
    int_high,
    event_time,
    max_iter = 100000,
    dgt = 3,
    seed) {
  
  iter <- 1
  while(abs(prop.cens.empirical - prop.cens.target) > tolerance && iter <= max_iter) {
    set.seed(seed)
    int_mid <- (int_low + int_high) / 2
    
    # Simulating event times
    # u <- runif(1000)
    # x1 <- rnorm(1000)
    # beta1 <- log(1.5)
    # lp1 <- beta1 * x1
    # scale <- 0.01
    # shape <- 2
    # event_time <- (-log(u)/ (scale * exp(lp1)))**(1/shape)
    
    # Simulating censoring
    cens <- rexp(n, rate = int_mid)
    
    # Observed times
    time <- pmin(event_time, cens)
    status <- as.numeric(event_time <= cens)
    prop.cens.empirical <- 1 - mean(status)
    
    # Print output
    # cat(sprintf("Iter %d: rate = %.4f, cens = %.4f\n", iter, int_mid, prop.cens.empirical))
    
    # Bisection
    if(prop.cens.empirical < prop.cens.target) {
      int_low <- int_mid
    } else {
      int_high <- int_mid
    }
    
    iter <- iter + 1
  }
  cat(sprintf("\nFinal rate = %.4f (target = %.2f, empirical = %.4f)\n", 
              int_mid, 
              prop.cens.target, 
              prop.cens.empirical))
  
  # Save results
  res_rate <- data.frame(
    "final_rate" = int_mid,
    "cens_target" = prop.cens.target,
    "cens_empirical" = prop.cens.empirical)
  
  res_rate <- round(res_rate, dgt)
  
  # Generate output
  # Censoring with the specified rate
  cens <- rexp(n, 
               rate = res_rate$final_rate)
  
   df <- data.frame(
     event_time = event_time,
     time = pmin(event_time, cens),
     status = as.numeric(event_time <= cens))

    res <- list(
      "res_rate" = res_rate,
      "df_output" = df)
    
    return(res)}
  





 
