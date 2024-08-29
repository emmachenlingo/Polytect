#' Plotting function for clustering results
#' 
#' This function takes results from polytect_clust and polytect_merge, or a data frame containing flurescence intensities and partition 
#' labels. It will output all combination of 2-color plots.
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_clust} and \code{polytect_merge} or any data frame containing the above information.
#' @return 2-color plots.
#' @examples
#' data(HR)
#' df_data<-polytect_clust(HR,4)
#' polytect_plot(df_data)
#' @export
polytect_plot<-function(df_data, cluster_num, cluster_selected=TRUE){
    col_num<-ncol(df_data)
    mat_select<-cluster_selection(cluster_num)
    
    plots<-list()
    k=0
    for (i in seq_len(col_num-2)){
        col_seq<-seq(i+1,col_num-1)
        for (j in col_seq){
            k=k+1
            df_data_tmp<-df_data[,c(i,j,col_num)]
            if(cluster_selected & col_num>3){
                if (col_num==4){
                    row_tmp<-which(mat_select[,-c(i,j)]==1)
                } else{
                    row_tmp<-apply(mat_select[,-c(i,j)],2,function(x) which(x==1))
                }
                row_selected<-unique(row_tmp)
                df_data_tmp<-df_data_tmp[!(df_data_tmp$cluster %in% row_selected),]
            }
            x_col <- colnames(df_data_tmp)[1]
            y_col <- colnames(df_data_tmp)[2]
            
            plots[[k]]<-ggplot(data=df_data_tmp, aes_string(x=x_col, y=y_col, colour = 'factor(cluster)'))+
                geom_point(size=0.9,show.legend = FALSE) +labs(x = paste("color",i), y=paste("color",j))+theme(text = element_text(size = 15),
                                                                                                               panel.grid.major = element_blank(),
                                                                                                               panel.grid.minor = element_blank(),
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(colour = "black"),
                                                                                                               plot.margin = margin(t = 20,  # Top margin
                                                                                                                                    r = 20,  # Right margin
                                                                                                                                    b = 20,  # Bottom margin
                                                                                                                                    l = 20))
            
        }
    }
    
    for (i in seq_len(k)) {
        plots[[i]] <- plots[[i]] +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  legend.position = "none",
                  plot.margin = unit(c(0.3, 0, 0, 0), "cm"))
    }
    
    p <- plot_grid(plotlist = plots, nrow = ceiling(k / 3), labels = LETTERS[1:k], label_size = 15, align = 'vh', hjust = 0)
    print(p)
}
