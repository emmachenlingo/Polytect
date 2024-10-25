#' Function for merging
#'
#' This function takes the clustering result as input. Users can first perform any clustering algorithm, then use this function. It
#' will return a data frame of fluorescence intensities and partition labels.
#' 
#' @importFrom stats cov setNames
#' @import utils
#' @import grDevices
#' @importFrom flowPeaks flowPeaks
#' @import ggplot2
#' @importFrom mvtnorm dmvnorm
#' @import dplyr
#' @import tidyverse
#' @importFrom cowplot plot_grid
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination setMBOControlInfill mbo
#' @importFrom DiceKriging km
#' @importFrom sn tr
#' @importFrom smoof makeSingleObjectiveFunction
#' @importFrom ParamHelpers makeParamSet generateDesign
#' @importFrom lhs randomLHS maximinLHS
#' @importFrom rgenoud genoud
#' @importFrom BiocManager install
#' @param data A matrix of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters.
#' @param base_clust A list that contains partition labels given by initial clustering.
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A data frame containing the original fluorescence intensity and the cluster labels.
#' @examples
#' data(HR)
#' dist_matrix <- dist(HR)
#' hc <- hclust(dist_matrix, method = "ward.D2")
#' hc_clusters <- cutree(hc, k = 6)
#' base_clust<-list()
#' base_clust$cluster<-hc_clusters
#' head(polytect_merge(HR, 4, base_clust))
#' @export
polytect_merge<-function(data,cluster_num,base_clust,lambdas=rep(2,64-log2(64)),coefs=rep(1,6)){
    data_scaled<-apply(data,2,function(x) (x-min(x))/(max(x)-min(x)))
    data_input<-as.matrix(data_scaled)
    
    g_clusternum<-unique(base_clust$cluster)
    
    df_data<-as.data.frame(cbind(data_input,cluster=base_clust$cluster))
    cluster_centers <- df_data %>%
        group_by(.data$cluster) %>%
        summarise(across(seq_len(ncol(df_data) - 1), mean, na.rm = TRUE))
    
    base_clust$mu<-as.matrix(cluster_centers)[,-1]
    
    
    result <- HMM_merge(data_input, cluster_num = cluster_num, base_clust = base_clust, eps = 10^(-10), max_iter = 1000, lambdas = lambdas[seq_len(cluster_num - log2(cluster_num))], coefs = coefs[seq_len(log2(cluster_num))])
    
    result_class<-apply(result[[1]],1,which.max)
    
    # Use the recode function from dplyr to update the 'group' column
    new_group <- recode(base_clust$cluster, !!!setNames(result_class, seq_along(g_clusternum)))
    df_data<-cbind(data,cluster=new_group)
    
    column_names <- colnames(df_data)
    
    # Rename the first n-1 columns
    new_column_names <- c(paste0("channel", seq_len(length(column_names) - 1)), column_names[length(column_names)])
    
    # Assign the new column names to the dataframe
    colnames(df_data) <- new_column_names
    
    return(as.data.frame(df_data))
}
