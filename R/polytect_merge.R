#' General function for merging
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
#' @param fp_par The parameters for flowPeaks. \code{fp_par}=c("default","manual","auto"). When "default" is chosen, the default parameters of
#' flowPeaks will be used. With "manual", you have to fill in \code{fp_optim}.
#' @param fp_optim The paramters for flowPeaks that users have to fill in manually when \code{fp_par} is set at "manual".
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#'
#' @return A data frame containing the original fluorescence intensity and the cluster labels.
#' @examples
#' data(HR)
#' polytect_merge(HR, 4)
#' @export
polytect_merge<-function(data,cluster_num,type="2color",fp_par="default",fp_optim=c(0.1,1,1.5),lambdas=rep(2,12),coefs=rep(1,4)){
  data_scaled<-apply(data,2,function(x) (x-min(x))/(max(x)-min(x)))
  data_input<-as.matrix(data_scaled)

  fp_tmp<-flowPeaks(data_input,tol=0.1,h0=1,h=1.5)

  if (fp_par=="default"){
    fp<-fp_tmp
  } else if (fp_par=='manual'){
    fp<-flowPeaks(data_input,tol=fp_optim[1],h0=fp_optim[2],h=fp_optim[3])
  } else if (fp_par=='auto'){
    hpo_result<-fp_search(data_input,cluster_num=cluster_num)
    fp<-flowPeaks(data_input,tol=hpo_result[1],h0=hpo_result[2],h=hpo_result[3])
  } else{
    return(print("The parameters of flowPeaks were specified wrong."))
  }
  g_clusternum<-unique(fp$peaks.cluster)
  g_clusternum_tmp<-unique(fp_tmp$peaks.cluster)

  if (length(g_clusternum_tmp)==1){
    df_data<-cbind(data,group=fp_tmp$peaks.cluster)
    return(df_data)
  }

  if(type=='2color'){
    result<-HMM_merge(data_input,cluster_num=4,fp,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:2],coefs=coefs[1:2])
  } else if(type=="2colorHO"){
    result<-HMM_merge_higher_order(data_input,cluster_num=8,fp,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:5],coefs=coefs[1:3])
  } else if(type=="3color"){
    result<-HMM_merge_3d(data_input,cluster_num=8,fp,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:5],coefs=coefs[1:3])
  } else if(type=="4color"){
    result<-HMM_merge_4d(data_input,cluster_num=16,fp,eps=10^(-10),max_iter=1000,lambdas=lambdas[1:12],coefs=coefs[1:4])
  } else {
    return(print("Warning: wrong type of data"))
  }
  result_class<-apply(result[[1]],1,which.max)

  # Use the recode function from dplyr to update the 'group' column
  new_group = recode(fp$peaks.cluster, !!!setNames(result_class, 1:length(g_clusternum)))
  df_data<-cbind(data,group=new_group)

  column_names <- colnames(df_data)

  # Rename the first n-1 columns
  new_column_names <- c(paste0("channel", 1:(length(column_names) - 1)), column_names[length(column_names)])

  # Assign the new column names to the dataframe
  colnames(df_data) <- new_column_names

  return(as.data.frame(df_data))
}
