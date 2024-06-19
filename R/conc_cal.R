#' concentration calculation function
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_merge} or any data frame containing the above information.
#' @param sampvol The sample volume in microliters (µL)
#' @param volmix The volume of the mixture
#' @param voltemp The volume of the template
#' @param type The assay design, including the number of channels and targets. \code{type}=c("2color",
#' "2colorHO","3color","4color"). "2color" is chosen when there are 2 colors and 2 targets. "2colorHO" means higher-order 2-color data (2 color and 3 targets). "3color" means
#' 3-color and 3-target. "4-color" is chosen when there are 4 colors.
#' @return a data frame of the summary of cluster centers, cluster sizes and cluster silhouette coefficients.
#' @examples
#' data(HR)
#' df_data<-polytect_merge(HR,4)
#' conc_cal(df_data)
#' @export
conc_cal<-function(df_data,sampvol=0.91,volmix=20,voltemp=20,type="2color"){
  summary_df<-polytect_summary(df_data)
  par_n<-sum(summary_df$group_size)
  if(type=='2color'){
    x1<-summary_df[summary_df$group==2,'group_size']+summary_df[summary_df$group==4,'group_size']
    x2<-summary_df[summary_df$group==2,'group_size']+summary_df[summary_df$group==4,'group_size']
    target1<-(1000/sampvol * (-log(1-x1/par_n)))*(volmix/voltemp)
    target2<-(1000/sampvol * (-log(1-x2/par_n)))*(volmix/voltemp)

    df_conc<-data.frame(target=c("1","2"),concentration=c(target1,target2))

  } else if (type=="2colorHO"|type=='3color'){
    x1<-summary_df[summary_df$group==2,'group_size']+summary_df[summary_df$group==5,'group_size']+summary_df[summary_df$group==6,'group_size']+summary_df[summary_df$group==8,'group_size']
    x2<-summary_df[summary_df$group==3,'group_size']+summary_df[summary_df$group==5,'group_size']+summary_df[summary_df$group==7,'group_size']+summary_df[summary_df$group==8,'group_size']
    x3<-summary_df[summary_df$group==4,'group_size']+summary_df[summary_df$group==6,'group_size']+summary_df[summary_df$group==7,'group_size']+summary_df[summary_df$group==8,'group_size']

    target1<-(1000/sampvol * (-log(1-x1/par_n)))*(volmix/voltemp)
    target2<-(1000/sampvol * (-log(1-x2/par_n)))*(volmix/voltemp)
    target3<-(1000/sampvol * (-log(1-x3/par_n)))*(volmix/voltemp)
    df_conc<-data.frame(target=c("1","2","3"),concentration=c(target1,target2,target3))
  } else if (type=='4color'){
    x1<-summary_df[summary_df$group==2,'group_size']+summary_df[summary_df$group==6,'group_size']+summary_df[summary_df$group==7,'group_size']+summary_df[summary_df$group==8,'group_size']+
      summary_df[summary_df$group==12,'group_size']+summary_df[summary_df$group==13,'group_size']+summary_df[summary_df$group==14,'group_size']+summary_df[summary_df$group==16,'group_size']
    x2<-summary_df[summary_df$group==3,'group_size']+summary_df[summary_df$group==6,'group_size']+summary_df[summary_df$group==9,'group_size']+summary_df[summary_df$group==10,'group_size']+
      summary_df[summary_df$group==12,'group_size']+summary_df[summary_df$group==13,'group_size']+summary_df[summary_df$group==15,'group_size']+summary_df[summary_df$group==16,'group_size']
    x3<-summary_df[summary_df$group==4,'group_size']+summary_df[summary_df$group==7,'group_size']+summary_df[summary_df$group==9,'group_size']+summary_df[summary_df$group==11,'group_size']+
      summary_df[summary_df$group==12,'group_size']+summary_df[summary_df$group==14,'group_size']+summary_df[summary_df$group==15,'group_size']+summary_df[summary_df$group==16,'group_size']
    x4<-summary_df[summary_df$group==5,'group_size']+summary_df[summary_df$group==8,'group_size']+summary_df[summary_df$group==10,'group_size']+summary_df[summary_df$group==11,'group_size']+
      summary_df[summary_df$group==13,'group_size']+summary_df[summary_df$group==14,'group_size']+summary_df[summary_df$group==15,'group_size']+summary_df[summary_df$group==16,'group_size']

    target1<-(1000/sampvol * (-log(1-x1/par_n)))*(volmix/voltemp)
    target2<-(1000/sampvol * (-log(1-x2/par_n)))*(volmix/voltemp)
    target3<-(1000/sampvol * (-log(1-x3/par_n)))*(volmix/voltemp)
    target4<-(1000/sampvol * (-log(1-x4/par_n)))*(volmix/voltemp)
    df_conc<-data.frame(target=c("1","2","3","4"),concentration=c(target1,target2,target3,target4))
  }

  return(df_conc)
}
