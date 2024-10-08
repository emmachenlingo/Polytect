% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_data.R
\docType{data}
\name{HIV}
\alias{HIV}
\title{HIV data}
\format{
A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.

\describe{
\item{channel1}{fluorescence intensities of color 1}
\item{channel2}{fluorescence intensities of color 2}
\item{channel3}{fluorescence intensities of color 3}
\item{channel4}{fluorescence intensities of color 4}
}
}
\source{
\url{https://www.biorxiv.org/content/10.1101/2023.08.18.553846v1}
}
\usage{
data(HIV)
}
\description{
A 4-color dPCR data of intact HIV-1 proviruses
}
\examples{
data(HIV)
head(HIV)
}
\keyword{datasets}
