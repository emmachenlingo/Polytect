#' concentration calculation function
#'
#' This function takes a data frame of fluorescence intensities and partition clusters as input. It can be results from polytect_clust or
#' polytect_merge. It will give the target concentration as output.
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_merge} or any data frame containing the above information.
#' @param cluster_num the expected number of clusters
#' @param sampvol The sample volume in microliters (µL)
#' @param volmix The volume of the mixture
#' @param voltemp The volume of the template
#' @return a data frame of target concentration.
#' @examples
#' data(HR)
#' df_data<-polytect_clust(HR,4)
#' conc_cal(df_data,4)
#' @export
conc_cal<-function(df_data,cluster_num,sampvol=0.91,volmix=20,voltemp=20){
    mat_coef<-cluster_selection(cluster_num)
    summary_df<-polytect_summary(df_data)
    par_n<-sum(summary_df$cluster_size)
    
    pos_clus<-apply(mat_coef,2,function(x) which(x==1))
    pos_pars<-apply(pos_clus,2,function(x) sum(summary_df[summary_df$cluster%in%x,'cluster_size']))
   
    targets<-(1000/sampvol * (-log(1-pos_pars/par_n)))*(volmix/voltemp)

    df_conc<-data.frame(target=as.character(seq_along(targets)),concentration=targets)
    
    return(df_conc)
}
