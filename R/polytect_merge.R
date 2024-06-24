#' Function for merging
#'
#' This function takes the clustering result as input. Users can first perform any clustering algorithm, then use this function. It
#' will return a data frame of fluorescence intensities and partition labels.
#' 
#' @import flowPeaks
#' @import ggplot2
#' @import mvtnorm
#' @import sn
#' @import dplyr
#' @import tidyverse
#' @import cowplot
#' @import mlrMBO
#' @import DiceKriging
#' @param data A matrix of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters.
#' @param type The assay design, including the number of channels and targets. \code{type}=c("2color",
#' "2colorHO","3color","4color"). "2color" is chosen when there are 2 colors and 2 targets. "2colorHO" means higher-order 2-color data (2 color and 3 targets). "3color" means
#' 3-color and 3-target. "4-color" is chosen when there are 4 colors.
#' @param base_clust A list that contains partition labels given by initial clustering.

#' @return A data frame containing the original fluorescence intensity and the cluster labels.
#' @examples
#' data(HR)
#' dist_matrix <- dist(HR)
#' hc <- hclust(dist_matrix, method = "ward.D2")
#' hc_clusters <- cutree(hc, k = 6)
#' polytect_merge(HR, 4, hc_clusters)
#' @export
polytect_merge<-function(data,cluster_num,base_clust,type="2color",lambdas=rep(2,12),coefs=rep(1,4)){
  data_scaled<-apply(data,2,function(x) (x-min(x))/(max(x)-min(x)))
  data_input<-as.matrix(data_scaled)
  
  g_clusternum<-unique(base_clust$cluster)

  df_data<-as.data.frame(cbind(data_input,cluster=base_clust$cluster))
  cluster_centers <- df_data %>%
    group_by(cluster) %>%
    summarise(across(1:(ncol(df_data)-1), mean, na.rm = TRUE))
  
  base_clust$mu<-as.matrix(cluster_centers)[,-1]
  

  if(type=='2color'){
    result<-HMM_merge(data_input,cluster_num=4,base_clust=base_clust,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:2],coefs=coefs[1:2])
  } else if(type=="2colorHO"){
    result<-HMM_merge_higher_order(data_input,cluster_num=8,base_clust=base_clust,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:5],coefs=coefs[1:3])
  } else if(type=="3color"){
    result<-HMM_merge_3d(data_input,cluster_num=8,base_clust=base_clust,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:5],coefs=coefs[1:3])
  } else if(type=="4color"){
    result<-HMM_merge_4d(data_input,cluster_num=16,base_clust=base_clust,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:12],coefs=coefs[1:4])
  } else {
    return(print("Warning: wrong type of data"))
  }
  result_class<-apply(result[[1]],1,which.max)
  
  # Use the recode function from dplyr to update the 'group' column
  new_group = recode(base_clust$cluster, !!!setNames(result_class, 1:length(g_clusternum)))
  df_data<-cbind(data,cluster=new_group)
  
  column_names <- colnames(df_data)
  
  # Rename the first n-1 columns
  new_column_names <- c(paste0("channel", 1:(length(column_names) - 1)), column_names[length(column_names)])
  
  # Assign the new column names to the dataframe
  colnames(df_data) <- new_column_names
  
  return(as.data.frame(df_data))
}
