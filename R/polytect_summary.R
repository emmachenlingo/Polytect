#' summary function
#'
#' @param df_data A data frame containing partition fluorescence intensities and corresponding cluster label. This can be the output
#' of \code{polytect_clust} or any data frame containing the above information.
#' @return a data frame of the summary of cluster centers, cluster sizes and cluster silhouette coefficients.
#' @examples
#' data(HR)
#' df_data<-polytect_clust(HR,4)
#' polytect_summary(df_data)
#' @export
polytect_summary<-function(df_data){
  # Identify the channel columns dynamically
  data_scaled<-apply(df_data[,-ncol(df_data)],2,function(x) (x-min(x))/(max(x)-min(x)))
  data_input<-as.matrix(data_scaled)
  df_data2<-data.frame(cbind(data_input,cluster=df_data$cluster))

  sil_coefs<-silhouette_coef(df_data2,df_data2$cluster)

  channel_columns <- df_data %>% select(starts_with("channel")) %>% colnames()

  # Group by 'group' and calculate the mean for each channel
  summary_df <- df_data %>%
    group_by(cluster) %>%
    summarise(across(all_of(channel_columns), mean, .names = "mean_{col}"),
              cluster_size = n())

  summary_df <- cbind(summary_df,silhouette_coef=sil_coefs[[2]][,2])

  return(summary_df)
}
