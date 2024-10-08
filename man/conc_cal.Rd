% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conc_cal.R
\name{conc_cal}
\alias{conc_cal}
\title{concentration calculation function}
\usage{
conc_cal(df_data, cluster_num, sampvol = 0.91, volmix = 20, voltemp = 20)
}
\arguments{
\item{df_data}{A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
of \code{polytect_merge} or any data frame containing the above information.}

\item{cluster_num}{the expected number of clusters}

\item{sampvol}{The sample volume in microliters (µL)}

\item{volmix}{The volume of the mixture}

\item{voltemp}{The volume of the template}
}
\value{
a data frame of target concentration.
}
\description{
This function takes a data frame of fluorescence intensities and partition clusters as input. It can be results from polytect_clust or
polytect_merge. It will give the target concentration as output.
}
\examples{
data(HR)
df_data<-polytect_clust(HR,4)
conc_cal(df_data,4)
}
