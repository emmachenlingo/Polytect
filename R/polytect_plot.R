#' plotting function
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_merge} or any data frame containing the above information.
#' @return 2-d plots.
#' @examples
#' data(HR)
#' df_data<-polytect_merge(HR,4)
#' polytect_plot(df_data)
#' @export
polytect_plot<-function(df_data){
  col_num<-ncol(df_data)
  plots<-list()
  k=0
  for (i in 1:(col_num-2)){
    for (j in (i+1):(col_num-1)){
      k=k+1
      df_data_tmp<-df_data[,c(i,j,col_num)]
      x_col <- colnames(df_data_tmp)[1]
      y_col <- colnames(df_data_tmp)[2]

      plots[[k]]<-ggplot(data=df_data_tmp, aes_string(x=x_col, y=y_col, colour = 'factor(group)'))+
        geom_point(size=0.9,show.legend = FALSE) +labs(x = paste("channel",i), y=paste("channel",j))+theme(text = element_text(size = 15),
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

  for (i in 1:k) {
    plots[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "none",
            plot.margin = unit(c(0.3, 0, 0, 0), "cm"))
  }

  p <- plot_grid(plotlist = plots, nrow = ceiling(k / 3), labels = LETTERS[1:k], label_size = 15, align = 'vh', hjust = 0)
  print(p)

}
