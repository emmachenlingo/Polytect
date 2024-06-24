# R/my_data.R
#' HR data
#'
#' A high-resolution 2-d dPCR data
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' good separation but some crosstalk; RPP30 genomic DNA assay
#' \describe{
#'   \item{channel1}{channel 1}
#'   \item{channel2}{channel 2}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/33992770/}
"HR"


#' MM data
#'
#' A multi-mode 2-d dPCR data
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' obvious multimodality; HIV gBlock sequences
#' \describe{
#'   \item{channel1}{channel 1}
#'   \item{channel2}{channel 2}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/37827643/}
"MM"


#' LR data
#'
#' A low-resolution 2-d dPCR data
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' barely separable on x-axis; development of genotyping assays for plants various; primers/probes and conditions for each target are evaluated.
#' \describe{
#'   \item{channel1}{channel 1}
#'   \item{channel2}{channel 2}
#' }
"LR"



