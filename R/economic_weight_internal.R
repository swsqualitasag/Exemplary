
#' @title
#' Internal Function using distribution thresholds.
#' @description
#' This function creates an economic value for a trait using the distribution thresholds of the traits records over the population. These thresholds are defined by price differences in a specific range of trait values. Price changes per kg carcass weight have been fitted to these thresholds. The difference of price changes between the initial distribution and the distribution with an increased mean divided by the differences in means is the economic value.
#' @param pvec_threshold   thresholds under continuous distribution given by pricing system
#' @param pvec_price       vector of prices for given classes
#' @param pn_mean          current population phenotypic mean
#' @param pn_sd            current population phenotypic standard deviation
#' @param pn_delta_mean    small change of current phenotypic population mean
#' @return                 economic value
#' @export compute_internal_economic_value
compute_internal_economic_value <- function(pvec_threshold,
                                            pvec_price,
                                            pn_mean,
                                            pn_sd,
                                            pn_delta_mean){
  # define new mean
  pn_mean_new <- pn_mean+pn_delta_mean
  # define frequencies for initial mean
  freq_pn_mean <- diff(c(0,pnorm(pvec_threshold,sd=pn_sd,mean=pn_mean, lower.tail = T)))
  # define frequencies for new mean
  freq_pn_mean_new <- diff(c(0,pnorm(pvec_threshold,sd=pn_sd,mean=pn_mean_new, lower.tail = T)))
  # compute economic weight
  ev_result <- (t(pvec_price)%*%freq_pn_mean_new-t(pvec_price)%*%freq_pn_mean)/pn_delta_mean
  return(ev_result)
}
