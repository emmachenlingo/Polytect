#' Main function for clustering
#' 
#' This is the main function for clustering. The function will start with flowPeaks, then merge the excess clusters. It will
#' return a data frame of fluorescence intensities and partition labels.
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
#' @importFrom smoof makeSingleObjectiveFunction
#' @importFrom sn tr
#' @importFrom ParamHelpers makeParamSet generateDesign
#' @importFrom lhs randomLHS maximinLHS
#' @importFrom rgenoud genoud
#' @importFrom BiocManager install
#' @param data A matrix of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters.
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
#' head(polytect_clust(HR, 4))
#' @export
polytect_clust<-function(data,cluster_num,fp_par="default",fp_optim=c(0.1,1,1.5),lambdas=rep(2,64-log2(64)),coefs=rep(1,6)){
    data_scaled<-apply(data,2,function(x) (x-min(x))/(max(x)-min(x)))
    data_input<-as.matrix(data_scaled)
    
    fp_tmp<-flowPeaks(data_input,tol=0.1,h0=1,h=1.5)
    
    if(length(unique(fp_tmp$kmeans.cluster)) < cluster_num){
        min_val <- apply(data_input, 2, min)
        dist_orig <- apply(data_input, 1, function(x) sqrt(sum((x-min_val)^2)))
        neg_ind <- which.min(dist_orig)
        repeated_row <- data_input[rep(neg_ind, 2*nrow(data_input)), ]
        data_input<-rbind(data_input,repeated_row)
    }
    
    if (fp_par=="default"){
        fp<-fp_tmp
    } else if (fp_par=='manual'){
        fp<-flowPeaks(data_input,tol=fp_optim[1],h0=fp_optim[2],h=fp_optim[3])
    } else if (fp_par=='auto'){
        hpo_result<-fp_search(data_input,cluster_num=cluster_num)
        fp<-flowPeaks(data_input,tol=hpo_result[1],h0=hpo_result[2],h=hpo_result[3])
    } else{
        stop("The parameters of flowPeaks were specified wrong.")
    }
    g_clusternum<-unique(fp$peaks.cluster)
    g_clusternum_tmp<-unique(fp_tmp$peaks.cluster)
    
    if (length(g_clusternum_tmp)==1){
        df_data<-cbind(data,cluster=fp_tmp$peaks.cluster)
        return(df_data)
    }
    fp_parse<-list()
    fp_parse$cluster<-fp$peaks.cluster[seq_len(nrow(data_scaled))]
    fp_parse$mu<-fp$peaks$mu
    
    data_input<-data_input[seq_len(nrow(data_scaled)),]
    result <- HMM_merge(data_input, cluster_num = cluster_num, base_clust = fp_parse, eps = 10^(-10), max_iter = 1000, lambdas = lambdas[seq_len(cluster_num - log2(cluster_num))], coefs = coefs[seq_len(log2(cluster_num))])
    result_class<-apply(result[[1]],1,which.max)
    # Use the recode function from dplyr to update the 'group' column
    new_group <- recode(fp_parse$cluster, !!!setNames(result_class, seq_along(g_clusternum)))
    df_data<-cbind(data,cluster=new_group)
    column_names <- colnames(df_data)
    
    # Rename the first n-1 columns
    new_column_names <- c(paste0("channel", seq_len(length(column_names) - 1)), column_names[length(column_names)])
    # Assign the new column names to the dataframe
    colnames(df_data) <- new_column_names
    
    return(as.data.frame(df_data))
}
