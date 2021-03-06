% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance.R
\name{compute_covariance_matrices}
\alias{compute_covariance_matrices}
\title{Covariance matrices between true and estimated breeding values.}
\usage{
compute_covariance_matrices(genetic_var_cov, residual_var_cov, offspring,
  proportion_calves, proportion_adults)
}
\arguments{
\item{genetic_var_cov}{The genetic variance-covariance matrix between the involved traits in the aggregate genotype.}

\item{residual_var_cov}{The genetic variance-covariance matrix between the involved traits in the aggregate genotype.}

\item{proportion_calves}{Proportion of calves in the offspring.}

\item{proportion_adults}{Proportion of adults in the offspring.}

\item{male_offspring}{Number of the assumed offspring of the male selection candidate when selection occurs.}

\item{female_offspring}{Number of the assumed offspring of the female selection candidate when selection occurs.}
}
\value{
Covariance Matrices between true and estimated breeding values
}
\description{
This function creates two covariance matrices between true and estimated breeding values of a male and a female selection candidate. These matrices are used to compute the genetic gain of an aggregate genotype when selection is based on multiple breeding values. The computation of the matrix is based on a possible scenario with known numbers of related animals, known genetic and residual variances and covariances. Here it is assumed, that there exist 6 breeding values which are present in the covariance matrix. This should result in two 6 x 6 matrices. Special is that it is also assumed, that the selection candidate has no records itself and that its relatives can have maximally have records of three of the six traits. This because they can not have both a record as calf and as adult. The variance-covariance matrices need the order calf, adult, calf, adult, calf, adult. In the result the matrices are bond together to one matrix, where the male covariance matrix is the matrix from row 1 to 6 and the female matrix is the matrix from row 7 to 12.
}
