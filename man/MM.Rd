% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_data.R
\docType{data}
\name{MM}
\alias{MM}
\title{MM data}
\format{
A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
obvious multimodality.
\describe{
\item{channel1}{fluorescence intensities of color 1}
\item{channel2}{fluorescence intensities of color 2}
}
}
\source{
\url{https://pubmed.ncbi.nlm.nih.gov/37827643/}
}
\usage{
data(MM)
}
\description{
A multi-mode 2-color dPCR data of HIV gBlock sequences
}
\examples{
data(MM)
head(MM)
}
\keyword{datasets}
