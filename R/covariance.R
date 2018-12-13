
#' @title
#' Covariance matrices between true and estimated breeding values.
#' @description
#' This function creates two covariance matrices between true and estimated breeding values of a male and a female selection candidate. These matrices are used to compute the genetic gain of an aggregate genotype when selection is based on multiple breeding values. The computation of the matrix is based on a possible scenario with known numbers of related animals, known genetic and residual variances and covariances. Here it is assumed, that there exist 6 breeding values which are present in the covariance matrix. This should result in two 6 x 6 matrices. Special is that it is also assumed, that the selection candidate has no records itself and that its relatives can have maximally have records of three of the six traits. This because they can not have both a record as calf and as adult. The variance-covariance matrices need the order calf, adult, calf, adult, calf, adult. In the result the matrices are bond together to one matrix, where the male covariance matrix is the matrix from row 1 to 6 and the female matrix is the matrix from row 7 to 12.
#' @param genetic_var_cov   The genetic variance-covariance matrix between the involved traits in the aggregate genotype.
#' @param residual_var_cov  The genetic variance-covariance matrix between the involved traits in the aggregate genotype.
#' @param male_offspring    Number of the assumed offspring of the male selection candidate when selection occurs.
#' @param female_offspring    Number of the assumed offspring of the female selection candidate when selection occurs.
#' @param proportion_calves Proportion of calves in the offspring.
#' @param proportion_adults Proportion of adults in the offspring.
#' @return                 Covariance Matrices between true and estimated breeding values
#' @export compute_covariance_matrices
compute_covariance_matrices <- function(genetic_var_cov,
                                        residual_var_cov,
                                        offspring,
                                        proportion_calves,
                                        proportion_adults){
  # define number of traits
  n <- as.numeric(nrow(genetic_var_cov))

  # Compute the relationsship matrix dependent on the number of offspring of the selection candidate. The selection candidate is the first animal. Male
  suppressPackageStartupMessages(library(pedigreemm))
  mnumb <- male_offspring+1
  male_A <- pedigree(sire = c(NA,rep(1,times = male_offspring)), dam =c(NA,rep(NA,times = male_offspring)), label = 1:mnumb)
  # Compute inverse of relationsship matrix.
  male_Ainv <- getAInv(male_A)

  # Compute the relationsship matrix dependent on the number of offspring of the selection candidate. The selection candidate is the first animal. Female
  suppressPackageStartupMessages(library(pedigreemm))
  mnumb <- female_offspring+1
  female_A <- pedigree(sire = c(NA,rep(NA,times = female_offspring)), dam =c(NA,rep(1,times = female_offspring)), label = 1:mnumb)
  # Compute inverse of relationsship matrix.
  female_Ainv <- getAInv(female_A)

  # Compute Relationship Matrix Z
  # numbers of adults and calves in scenario
  male_number_adults <- floor(proportion_adults*male_offspring)[[1]]
  male_number_calves <- male_offspring-male_number_adults
  female_number_adults <- floor(proportion_adults*female_offspring)[[1]]
  female_number_calves <- female_offspring-female_number_adults
  # Kronecker
  kronecker_calves <- cbind(c(1,0,0),c(0,0,0),c(0,1,0),c(0,0,0),c(0,0,1),c(0,0,0))
  kronecker_adults <- cbind(c(0,0,0),c(1,0,0),c(0,0,0),c(0,1,0),c(0,0,0),c(0,0,1))
  # male Z
  male_calves <- diag(1,male_number_calves)%x%kronecker_calves
  male_calves <- cbind(male_calves,matrix(0, nrow=nrow(male_calves),ncol=male_number_adults*n))
  male_adults <- diag(male_number_adults)%x%kronecker_adults
  male_adults <- cbind(matrix(0, nrow=nrow(male_adults),ncol=male_number_calves*n),male_adults)
  male_Z <- rbind(male_calves, male_adults)
  male_Z<-cbind(matrix(0, nrow=male_offspring*(n/2),ncol=n),male_Z)
  # female Z
  female_calves <- diag(1,female_number_calves)%x%kronecker_calves
  female_calves <- cbind(female_calves,matrix(0, nrow=nrow(female_calves),ncol=female_number_adults*n))
  female_adults <- diag(female_number_adults)%x%kronecker_adults
  female_adults <- cbind(matrix(0, nrow=nrow(female_adults),ncol=female_number_calves*n),female_adults)
  female_Z <- rbind(female_calves, female_adults)
  female_Z<-cbind(matrix(0, nrow=female_offspring*(n/2),ncol=n),female_Z)

  # Genetic Variance Covariance Matrix G
  G <- genetic_var_cov

  # Residual Variance Covariance Matrix R
  residual_var_cov_calves <- residual_var_cov[c(1,3,5),c(1,3,5)]
  residual_var_cov_adults <- residual_var_cov[c(2,4,6),c(2,4,6)]
  # male
  male_calves_kronecker_residual <- diag(male_number_calves) %x% residual_var_cov_calves
  male_adults_kronecker_residual <- diag(male_number_adults) %x% residual_var_cov_adults
  male_calves_kronecker_residual_extended <- cbind(male_calves_kronecker_residual, matrix(0,nrow=nrow(male_calves_kronecker_residual),ncol=ncol(male_adults_kronecker_residual)))
  male_adults_kronecker_residual_extended <- cbind( matrix(0,nrow=nrow(male_adults_kronecker_residual),ncol=ncol(male_calves_kronecker_residual)), male_adults_kronecker_residual)
  male_R <- rbind(male_calves_kronecker_residual_extended,male_adults_kronecker_residual_extended)
  # female
  female_calves_kronecker_residual <- diag(female_number_calves) %x% residual_var_cov_calves
  female_adults_kronecker_residual <- diag(female_number_adults) %x% residual_var_cov_adults
  female_calves_kronecker_residual_extended <- cbind(female_calves_kronecker_residual, matrix(0,nrow=nrow(female_calves_kronecker_residual),ncol=ncol(female_adults_kronecker_residual)))
  female_adults_kronecker_residual_extended <- cbind( matrix(0,nrow=nrow(female_adults_kronecker_residual),ncol=ncol(female_calves_kronecker_residual)), female_adults_kronecker_residual)
  female_R <- rbind(female_calves_kronecker_residual_extended,female_adults_kronecker_residual_extended)

  # Prediction error variance matrix
  # male
  male_PEV <- solve(t(male_Z)%*%solve(male_R)%*%male_Z+male_Ainv%x%solve(G))[1:n,1:n]
  # female
  female_PEV <- solve(t(female_Z)%*%solve(female_R)%*%female_Z+female_Ainv%x%solve(G))[1:n,1:n]

  # Covariance Matrix
  # male
  male_C <- G-male_PEV
  # female
  female_C <- G-female_PEV
  # together
  C <- rbind(male_C,female_C)

  # Return
  return(C)
}
