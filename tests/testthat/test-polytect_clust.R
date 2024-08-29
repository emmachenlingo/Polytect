## test the main function polytect_clust
library(testthat)
library(Polytect)

test_that("polytect_clust works as expected", {
    data("CNV6plex")
    df_data<-polytect_clust(CNV6plex,cluster_num=64,fp_par="auto",fp_optim=c(0.1,1,1.5),lambdas=rep(2,64-log2(64)),coefs=rep(1,6))
    expect_equal(length(unique(df_data$cluster)), 64)
})
