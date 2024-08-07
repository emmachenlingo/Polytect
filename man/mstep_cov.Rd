% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_functions.R
\name{mstep_cov}
\alias{mstep_cov}
\title{Internal Function 9}
\usage{
mstep_cov(cluster_num, dim_data, g_clusternum, zi, mg, covg, mug, muh)
}
\arguments{
\item{cluster_num}{The expected maximum number of clusters}

\item{dim_data}{the dimension of the dataset}

\item{g_clusternum}{cluster labels from base clustering}

\item{zi}{the expected log-likelihood found on the E step}

\item{mg}{cluster sizes of base clustering result}

\item{covg}{the covariance matrix of clusters at level l+1}

\item{mug}{the matrix of cluster centers at level l+1}

\item{muh}{the matrix of cluster centers at level l}
}
\value{
covh the covariance matrix of clusters at level l in the EM algorithm
}
\description{
This function calculates mu in M-step of EM algorithm
}
\keyword{internal}
