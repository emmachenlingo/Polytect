#' plotting function
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_merge} or any data frame containing the above information.
#' @return plot of silhouette coefficients for each cluster.
#' @examples
#' data(HR)
#' df_data<-polytect_merge(HR,4)
#' sil_plot(df_data)
#' @export
sil_plot<-function(df_data){
  data_scaled<-apply(df_data[,-ncol(df_data)],2,function(x) (x-min(x))/(max(x)-min(x)))
  data_input<-as.matrix(data_scaled)
  df_data2<-data.frame(cbind(data_input,cluster=df_data$cluster))

  sil_coefs<-silhouette_coef(df_data2,df_data2$cluster)


  # Plot the silhouette coefficients using ggplot2
  ggplot(sil_coefs[[1]], aes(x = cluster, y = width, group=factor(cluster),fill = factor(cluster))) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "silhouette plot", x = "cluster", y = "silhouette width", fill = "cluster") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)))

}
