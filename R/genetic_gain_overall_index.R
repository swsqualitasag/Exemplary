
#' @title
#' Genetic Gain in Sfr. / year for index.
#' @description
#' This function computes the genetic gain in Sfr. / year of the selection on the index. In this index all traits are represented in the aggregate genotype.
#' @param genetic_var_cov   The genetic variance-covariance matrix between the involved traits in the aggregate genotype.
#' @param residual_var_cov  The genetic variance-covariance matrix between the involved traits in the aggregate genotype.
#' @param male_offspring    Number of the assumed offspring of the male selection candidate when selection occurs.
#' @param female_offspring    Number of the assumed offspring of the female selection candidate when selection occurs.
#' @param proportion_calves Proportion of calves in the offspring.
#' @param proportion_adults Proportion of adults in the offspring.
#' @param male_proportionselected Proportion selected of the offspring of a male selection candidate
#' @param female_proportionselected Proportion selected of the offspring of a female selection candidate
#' @param male_generationintervall Generation intervall of the male selection path.
#' @param female_generationintervall Generation intervall of the female selection path.
#' @param economic_weights A vector of economic weights with lenght of number of traits.
#' @return                 Genetic Gain in Sfr. per year
#' @export compute_genetic_gain_overall_index
compute_genetic_gain_overall_index <- function(genetic_var_cov,
                                               residual_var_cov,
                                               offspring,
                                               proportion_calves,
                                               proportion_adults,
                                               male_proportionselected,
                                               female_proportionselected,
                                               male_generationintervall,
                                               female_generationintervall,
                                               economic_weights){
  # define vector of economic weights
  a <- economic_weights

  # Intesity of selection
  male_i <- dnorm(qnorm(1-male_proportionselected))/male_proportionselected
  female_i <- dnorm(qnorm(1-female_proportionselected))/female_proportionselected

  # Asymptotic covariance matrices

  Ca <- compute_asymptotic_covariance_matrices(genetic_var_cov = genetic_var_cov, residual_var_cov = residual_var_cov, offspring = offspring, proportion_calves = proportion_calves, proportion_adults = proportion_adults,male_proportionselected = male_proportionselected,female_proportionselected = female_proportionselected)

  # Number of traits
  n <- as.numeric(ncol(Ca))

  # male
  male_Ca <- Ca[1:n,]
  # female
  female_Ca <- Ca[(1+n):(2*n),]

  # Genetic Gain
  Qoverall_index <- (male_i*sqrt(as.numeric(t(a)%*%male_Ca%*%a))+female_i*sqrt(as.numeric(t(a)%*%female_Ca%*%a)))/(male_generationintervall+female_generationintervall)

  return(Qoverall_index)
}
