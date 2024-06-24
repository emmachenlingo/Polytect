#' concentration calculation function
#'
#' This function takes a data frame of fluorescence intensities and partition clusters as input. It can be results from polytect_clust or
#' polytect_merge. It will give the target concentration as output.
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_merge} or any data frame containing the above information.
#' @param sampvol The sample volume in microliters (µL)
#' @param volmix The volume of the mixture
#' @param voltemp The volume of the template
#' @param type The assay design, including the number of channels and targets. \code{type}=c("2color",
#' "2colorHO","3color","4color"). "2color" is chosen when there are 2 colors and 2 targets. "2colorHO" means higher-order 2-color data (2 color and 3 targets). "3color" means
#' 3-color and 3-target. "4-color" is chosen when there are 4 colors.
#' @return a data frame of target concentration.
#' @examples
#' data(HR)
#' df_data<-polytect_clust(HR,4)
#' conc_cal(df_data)
#' @export
conc_cal<-function(df_data,sampvol=0.91,volmix=20,voltemp=20,type="2color"){
  summary_df<-polytect_summary(df_data)
  par_n<-sum(summary_df$cluster_size)
  if(type=='2color'){
    x1<-summary_df[summary_df$cluster==2,'cluster_size']+summary_df[summary_df$cluster==4,'cluster_size']
    x2<-summary_df[summary_df$cluster==2,'cluster_size']+summary_df[summary_df$cluster==4,'cluster_size']
    target1<-(1000/sampvol * (-log(1-x1/par_n)))*(volmix/voltemp)
    target2<-(1000/sampvol * (-log(1-x2/par_n)))*(volmix/voltemp)

    df_conc<-data.frame(target=c("1","2"),concentration=c(target1,target2))

  } else if (type=="2colorHO"|type=='3color'){
    x1<-summary_df[summary_df$cluster==2,'cluster_size']+summary_df[summary_df$cluster==5,'cluster_size']+summary_df[summary_df$cluster==6,'cluster_size']+summary_df[summary_df$cluster==8,'cluster_size']
    x2<-summary_df[summary_df$cluster==3,'cluster_size']+summary_df[summary_df$cluster==5,'cluster_size']+summary_df[summary_df$cluster==7,'cluster_size']+summary_df[summary_df$cluster==8,'cluster_size']
    x3<-summary_df[summary_df$cluster==4,'cluster_size']+summary_df[summary_df$cluster==6,'cluster_size']+summary_df[summary_df$cluster==7,'cluster_size']+summary_df[summary_df$cluster==8,'cluster_size']

    target1<-(1000/sampvol * (-log(1-x1/par_n)))*(volmix/voltemp)
    target2<-(1000/sampvol * (-log(1-x2/par_n)))*(volmix/voltemp)
    target3<-(1000/sampvol * (-log(1-x3/par_n)))*(volmix/voltemp)
    df_conc<-data.frame(target=c("1","2","3"),concentration=c(target1,target2,target3))
  } else if (type=='4color'){
    x1<-summary_df[summary_df$cluster==2,'cluster_size']+summary_df[summary_df$cluster==6,'cluster_size']+summary_df[summary_df$cluster==7,'cluster_size']+summary_df[summary_df$cluster==8,'cluster_size']+
      summary_df[summary_df$cluster==12,'cluster_size']+summary_df[summary_df$cluster==13,'cluster_size']+summary_df[summary_df$cluster==14,'cluster_size']+summary_df[summary_df$cluster==16,'cluster_size']
    x2<-summary_df[summary_df$cluster==3,'cluster_size']+summary_df[summary_df$cluster==6,'cluster_size']+summary_df[summary_df$cluster==9,'cluster_size']+summary_df[summary_df$cluster==10,'cluster_size']+
      summary_df[summary_df$cluster==12,'cluster_size']+summary_df[summary_df$cluster==13,'cluster_size']+summary_df[summary_df$cluster==15,'cluster_size']+summary_df[summary_df$cluster==16,'cluster_size']
    x3<-summary_df[summary_df$cluster==4,'cluster_size']+summary_df[summary_df$cluster==7,'cluster_size']+summary_df[summary_df$cluster==9,'cluster_size']+summary_df[summary_df$cluster==11,'cluster_size']+
      summary_df[summary_df$cluster==12,'cluster_size']+summary_df[summary_df$cluster==14,'cluster_size']+summary_df[summary_df$cluster==15,'cluster_size']+summary_df[summary_df$cluster==16,'cluster_size']
    x4<-summary_df[summary_df$cluster==5,'cluster_size']+summary_df[summary_df$cluster==8,'cluster_size']+summary_df[summary_df$cluster==10,'cluster_size']+summary_df[summary_df$cluster==11,'cluster_size']+
      summary_df[summary_df$cluster==13,'cluster_size']+summary_df[summary_df$cluster==14,'cluster_size']+summary_df[summary_df$cluster==15,'cluster_size']+summary_df[summary_df$cluster==16,'cluster_size']

    target1<-(1000/sampvol * (-log(1-x1/par_n)))*(volmix/voltemp)
    target2<-(1000/sampvol * (-log(1-x2/par_n)))*(volmix/voltemp)
    target3<-(1000/sampvol * (-log(1-x3/par_n)))*(volmix/voltemp)
    target4<-(1000/sampvol * (-log(1-x4/par_n)))*(volmix/voltemp)
    df_conc<-data.frame(target=c("1","2","3","4"),concentration=c(target1,target2,target3,target4))
  }

  return(df_conc)
}
