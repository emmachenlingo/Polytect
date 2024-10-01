# R/my_data.R
#' HR data
#'
#' A high-resolution 2-color dPCR data of RPP30 genomic DNA assay
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' good separation but some crosstalk.
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#' }
#' @usage data(HR)
#' @examples
#' data(HR)
#' head(HR)
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/33992770/}
"HR"


#' MM data
#'
#' A multi-mode 2-color dPCR data of HIV gBlock sequences
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' obvious multimodality.
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#' }
#' @usage data(MM)
#' @examples
#' data(MM)
#' head(MM)
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/37827643/}
"MM"


#' LR data
#'
#' A low-resolution 2-color dPCR data of development of genotyping assays for plants various
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' barely separable on x-axis.
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#' }
#' @usage data(LR)
#' @examples
#' data(LR)
#' head(LR)
"LR"


#' CA data
#'
#' 2-color competitive assay of competition BRAF V600E assay with 1% mutant
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' data is not orthogonal.
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#' }
#' @usage data(CA)
#' @examples
#' data(CA)
#' head(CA)
"CA"


#' BPV data
#'
#' A 3-color dPCR data of bovine papilloma virus assay
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' 
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#'   \item{channel3}{fluorescence intensities of color 3}
#' }
#' @usage data(BPV)
#' @examples
#' data(BPV)
#' head(BPV)
"BPV"


#' HIV data
#'
#' A 4-color dPCR data of intact HIV-1 proviruses
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' 
#' \describe{
#'    \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#'   \item{channel3}{fluorescence intensities of color 3}
#'   \item{channel4}{fluorescence intensities of color 4}
#' }
#' @usage data(HIV)
#' @examples
#' data(HIV)
#' head(HIV)
#' @source \url{https://www.biorxiv.org/content/10.1101/2023.08.18.553846v1}
"HIV"


#' CNV 5-plex data
#'
#' CNV 5-plex universal probes
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' 
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#'   \item{channel3}{fluorescence intensities of color 3}
#'   \item{channel4}{fluorescence intensities of color 4}
#'   \item{channel5}{fluorescence intensities of color 5}

#' }
#' @usage data(CNV5plex)
#' @examples
#' data(CNV5plex)
#' head(CNV5plex)
"CNV5plex"


#' CNV 6-plex data
#'
#' CNV 6-plex universal probes
#'
#' @format A data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' 
#' \describe{
#'   \item{channel1}{fluorescence intensities of color 1}
#'   \item{channel2}{fluorescence intensities of color 2}
#'   \item{channel3}{fluorescence intensities of color 3}
#'   \item{channel4}{fluorescence intensities of color 4}
#'   \item{channel5}{fluorescence intensities of color 5}
#'   \item{channel6}{fluorescence intensities of color 6}

#' }
#' @usage data(CNV6plex)
#' @examples
#' data(CNV6plex)
#' head(CNV6plex)
"CNV6plex"



