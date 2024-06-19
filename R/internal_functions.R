#' Internal Function 1
#'
#' This function outputs silhouette coefficients.
#'
#' @param data A data frame containing standardized partition fluorescence intensities and corresponding cluster label.
#' @param clustering cluster labels
#' @return A list of silhouette coefficients for each partition and the mean silhouette coefficients for each cluster.
#' @keywords internal
silhouette_coef<-function(data,clustering,plot=FALSE,plot_name='orig',sim='orig'){

  # si <- silhouette(clustering, dist(data, "euclidean"))
  ndim<-ncol(data)
  si <- approxSilhouette(data[,-ndim], clustering)

  sil_tab<-as.data.frame(si) %>% group_by(cluster) %>% summarise(mean_sil=round(mean(width),2))

  # find the mean of each cluster

  clust_pos<- data %>% group_by(group) %>% summarise(mean_x=mean(channel1),mean_y=mean(channel2))
  colnames(clust_pos)[1]<-'cluster'
  df_merged <- merge(x=sil_tab,y=clust_pos, by="cluster", all.x=TRUE)

  if(plot){
    png(file=paste0(plot_name,sim,'.png'),width = 500,height = 300)
    p<-ggplot(data=data, aes(channel1, channel2, colour = factor(group)))+
      geom_point(size=0.7,show.legend = FALSE) +labs(x = "Green Channel",y='Red Channel')+theme(panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),
                                                                                                panel.background = element_blank(),
                                                                                                axis.line = element_line(colour = "black"),
                                                                                                plot.margin = margin(t = 20,  # Top margin
                                                                                                                     r = 20,  # Right margin
                                                                                                                     b = 20,  # Bottom margin
                                                                                                                     l = 20)) + annotate("text", x=df_merged$mean_x, y=df_merged$mean_y, label= df_merged$mean_sil)

    print(p)
    dev.off()
  }

  return(list(si,df_merged))

}

#' Internal Function 2
#'
#' This function outputs silhouette coefficients.
#'
#' @param data A matrix of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param clustering cluster labels
#' @return A data frame of silhouette coefficients for each partition.
#' @keywords internal
approxSilhouette <- function(x, clusters) {
  x <- as.matrix(x)
  uclust <- sort(unique(clusters))
  averaged <- list(length(uclust))
  clust.var <- numeric(length(uclust))

  for (i in seq_along(uclust)) {
    current <- uclust[i]==clusters
    xcurrent <- x[current,,drop=FALSE]
    centroid <- colMeans(xcurrent)
    averaged[[i]] <- centroid
    clust.var[i] <- sum(colMeans(sweep(xcurrent, 2, centroid)^2))
  }

  self.dist <- other.dist <- rep(Inf, nrow(x))
  other.clust <- integer(nrow(x))
  tx <- t(x)

  for (i in seq_along(uclust)) {
    D <- sqrt(colSums((tx - averaged[[i]])^2) + clust.var[i])

    is.self <- uclust[i]==clusters
    self.dist[is.self] <- D[is.self]

    is.other <- !is.self
    other.D <- D[is.other]
    better <- other.D < other.dist[is.other]
    other.dist[is.other][better] <- other.D[better]
    other.clust[is.other][better] <- i
  }

  result<-data.frame(
    cluster=clusters,
    other=uclust[other.clust],
    width=(other.dist - self.dist)/pmax(other.dist, self.dist),
    row.names=rownames(x)
  )

  return(result)
}


#' Internal Function 3
#'
#' This function optimizes parameters of flowPeaks
#' @param data A matrix of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @return A vector containing the optimal parameters found by the algorithm
#' @keywords internal
fp_search<-function(data,cluster_num=16){

  # simple 2d objective function
  obj.fun.sil = makeSingleObjectiveFunction(
    fn = function(pars) {
      fp<-flowPeaks(data,tol=pars[1],h0=pars[2],h=pars[3])
      deviation<-cluster_num+2-length(unique(fp$peaks.cluster))
      tryCatch( {
        fp_sil_optimval<-deviation^2
      },error = function(e) {
        fp_sil_optimval<<-10
      })

      return(fp_sil_optimval)
    },
    par.set = makeParamSet(
      makeNumericParam("x1", lower = 0, upper = 1),
      makeNumericParam("x2", lower = 0.1, upper = 5),
      makeNumericParam("x3", lower = 0.1, upper = 5)
    )
  )

  # create base control object
  ctrl = makeMBOControl()
  # do three MBO iterations
  ctrl = setMBOControlTermination(ctrl, iters = 4L)
  # use 500 points in the focus search (should be sufficient for 2d)
  ctrl = setMBOControlInfill(ctrl, opt.focussearch.points = 500)
  # create initial design
  des = generateDesign(n = 15L, getParamSet(obj.fun.sil), fun = lhs::maximinLHS)
  # start mbo
  res = mbo(obj.fun.sil, design = des, control = ctrl)

  bo_sil_result<-unlist(res$x)

  return(bo_sil_result)
}


#' Internal Function 4
#'
#' This function merges 2-d dPCR data
#' @param data A matrix or data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @param eps the convergence threshold
#' @param max_iter maximum number of iterations
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A list of membership probability, cluster center, merging probability
#' @keywords internal
HMM_merge<-function(data,cluster_num,fp,eps=10^(-10),max_iter=1000,lambdas=rep(2,2),coefs=rep(1,2)) {
  change_ests=NULL
  data<-as.matrix(data)
  dim_data<-ncol(data)

  pih = rep(1/cluster_num,cluster_num)
  g_clusternum<-unique(fp$peaks.cluster)
  mg<-table(fp$peaks.cluster)
  mug<-fp$peaks$mu
  covg <- array(0,dim=c(dim_data,dim_data,length(g_clusternum)))
  for (i in 1:length(g_clusternum)){
    tryCatch( {
      covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    }, error = function(e) {
      covg[,,g_clusternum[i]] <<- diag(0.001,nrow=dim_data)
    })
  }

  #initialization
  ## find the negative population
  min_val<-apply(data,2,min)
  dist_orig<-apply(mug,1,function(x) {sqrt(sum((x-min_val)^2))})
  neg_assum<-mug[which.min(dist_orig),]

  muh<-matrix(0,nrow=cluster_num,ncol=dim_data)
  muh[1,]<-neg_assum
  muh[2,]<-c(0.75*(min(mug[,1])+max(mug[,1])),neg_assum[2])
  muh[3,]<-c(neg_assum[1],0.75*(min(mug[,2])+max(mug[,2])))
  muh[4,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+(1-coefs[1]-coefs[2])*muh[1,]
  covh <- array(cov(mug), dim = c(dim_data, dim_data, cluster_num))

  zi<-matrix(0,nrow=length(g_clusternum), ncol=cluster_num)

  ## start the EM algorithm

  for (j in 1:max_iter){
    if (j>=max_iter){
      print("Warning: the algorithm fails to converge")
    }
    ## E step:
    for (g in 1:length(g_clusternum)){
      for (k in 1:cluster_num){
        inv_zi_tmp<-0
        for (l in 1:cluster_num){
          tmp=exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
          inv_zi_tmp=inv_zi_tmp+exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
        }
        zi[g,k]=1/inv_zi_tmp
      }
    }


    ## M step:
    pih = apply(zi,2,sum)/length(g_clusternum)
    pih[which(pih==0)]<-(10^(-10))
    ## for the mus
    mu1_nom<-c(0,0)
    mu1_denom<-c(0,0)
    mu2_nom<-c(0,0)
    mu2_denom<-c(0,0)
    mu3_nom<-c(0,0)
    mu3_denom<-c(0,0)
    mu4_nom<-c(0,0)
    mu4_denom<-c(0,0)
    for (g in 1:length(g_clusternum)){
      mu1_nom<-mu1_nom+zi[g,1]*mg[g]*mug[g,]%*%solve(covh[,,1])
      mu1_denom<-mu1_denom+zi[g,1]*mg[g]*solve(covh[,,1])

      mu2_nom<-mu2_nom+zi[g,2]*mg[g]*mug[g,]%*%solve(covh[,,2])
      mu2_denom<-mu2_denom+zi[g,2]*mg[g]*solve(covh[,,2])

      mu3_nom<-mu3_nom+zi[g,3]*mg[g]*mug[g,]%*%solve(covh[,,3])
      mu3_denom<-mu3_denom+zi[g,3]*mg[g]*solve(covh[,,3])

      mu4_nom<-mu4_nom+zi[g,4]*mg[g]*mug[g,]%*%solve(covh[,,4])
      mu4_denom<-mu4_denom+zi[g,4]*mg[g]*solve(covh[,,4])

    }
    mu1_nom<-mu1_nom+2*lambdas[1]*neg_assum+2*lambdas[2]*(coefs[1]+coefs[2]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-muh[4,])
    mu1_denom<-mu1_denom+2*lambdas[1]*diag(1,nrow=2)+2*lambdas[2]*(coefs[1]+coefs[2]-1)^2*diag(1,nrow=dim_data)

    mu2_nom<-mu2_nom+2*lambdas[2]*coefs[1]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[4,]-coefs[2]*muh[3,])
    mu2_denom<-mu2_denom+2*lambdas[2]*(coefs[1])^2*diag(1,nrow=dim_data)

    mu3_nom<-mu3_nom+2*lambdas[2]*coefs[2]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[4,]-coefs[1]*muh[2,])
    mu3_denom<-mu3_denom+2*lambdas[2]*(coefs[2])^2*diag(1,nrow=dim_data)

    mu4_nom<-mu4_nom+2*lambdas[2]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-(coefs[1]+coefs[2]-1)*muh[1,])
    mu4_denom<-mu4_denom+2*lambdas[2]*diag(1,nrow=dim_data)

    muh[1,]=mu1_nom%*%solve(mu1_denom)
    muh[2,]=mu2_nom%*%solve(mu2_denom)
    muh[3,]=mu3_nom%*%solve(mu3_denom)
    muh[4,]=mu4_nom%*%solve(mu4_denom)

    ## for Sigmas
    sigma_nom<-matrix(rep(0,dim_data^2),nrow=dim_data)
    sigma_denom<-0
    for (k in 1:cluster_num) {
      for (g in 1:length(g_clusternum)){
        sigma_nom<-sigma_nom+zi[g,k]*mg[g]*covg[,,g]+zi[g,k]*mg[g]*(mug[g,]-muh[k,])%*%(t(mug[g,]-muh[k,]))
        sigma_denom<-sigma_denom+zi[g,k]*mg[g]
      }
      covh[,, k] <- sigma_nom/sigma_denom
    }

    # change the likelihood
    change_ests<-c(change_ests,sum(apply(zi,1,function(x) log(sum(x)))))
    ## stopping criteria
    if(j>1){

      if(sum(abs(muh-old_mu))<eps){break}
    }

    old_mu <- muh


  }

  return(list(zi,muh,pih))
}


#' Internal Function 5
#'
#' This function merges 3-d dPCR data
#' @param data A matrix or data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @param eps the convergence threshold
#' @param max_iter maximum number of iterations
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A list of membership probability, cluster center, merging probability
#' @keywords internal
HMM_merge_3d<-function(data,cluster_num,fp,eps=10^(-10),max_iter=1000,lambdas=rep(2,5),coefs=rep(1,3)) {
  change_ests=NULL
  data<-as.matrix(data)
  dim_data<-ncol(data)

  pih = rep(1/cluster_num,cluster_num)
  g_clusternum<-unique(fp$peaks.cluster)
  mg<-table(fp$peaks.cluster)
  mug<-fp$peaks$mu
  covg <- array(0,dim=c(dim_data,dim_data,length(g_clusternum)))
  for (i in 1:length(g_clusternum)){
    # tryCatch( {
    #   covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    #   print(covg[,,g_clusternum[i]])
    # }, error = function(e) {
    #   covg[,,g_clusternum[i]] <<- covg[,,g_clusternum[i-1]]
    #   print(covg[,,g_clusternum[i-1]])
    # })
    tryCatch( {
      covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    }, error = function(e) {
      covg[,,g_clusternum[i]] <<- diag(0.001,nrow=dim_data)
    })
  }

  #initialization
  ## find the negative population
  min_val<-apply(data,2,min)
  dist_orig<-apply(mug,1,function(x) {sqrt(sum((x-min_val)^2))})
  neg_assum<-mug[which.min(dist_orig),]

  muh<-matrix(0,nrow=cluster_num,ncol=dim_data)
  muh[1,]<-neg_assum
  muh[2,]<-c(0.75*(min(mug[,1])+max(mug[,1])),neg_assum[2],neg_assum[3])
  muh[3,]<-c(neg_assum[1],0.75*(min(mug[,2])+max(mug[,2])),neg_assum[3])
  muh[4,]<-c(neg_assum[1],neg_assum[2],0.75*(min(mug[,3])+max(mug[,3])))
  muh[5,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+(1-coefs[1]-coefs[2])*muh[1,]
  muh[6,]<-coefs[1]*muh[2,]+coefs[3]*muh[4,]+(1-coefs[1]-coefs[3])*muh[1,]
  muh[7,]<-coefs[2]*muh[3,]+coefs[3]*muh[4,]+(1-coefs[2]-coefs[3])*muh[1,]
  muh[8,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]+(1-coefs[1]-coefs[2]-coefs[3])*muh[1,]

  covh <- array(cov(mug), dim = c(dim_data, dim_data, cluster_num))

  zi<-matrix(0,nrow=length(g_clusternum), ncol=cluster_num)

  ## start the EM algorithm

  for (j in 1:max_iter){
    if (j>=max_iter){
      print("Warning: the algorithm fails to converge")
    }
    ## E step:
    for (g in 1:length(g_clusternum)){
      for (k in 1:cluster_num){
        inv_zi_tmp<-0
        for (l in 1:cluster_num){
          tmp=exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
          inv_zi_tmp=inv_zi_tmp+exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
        }
        zi[g,k]=1/inv_zi_tmp
      }
    }


    ## M step:
    pih = apply(zi,2,sum)/length(g_clusternum)
    pih[which(pih==0)]<-(10^(-10))
    ## for the mus
    mu1_nom<-c(0,0,0)
    mu1_denom<-c(0,0,0)
    mu2_nom<-c(0,0,0)
    mu2_denom<-c(0,0,0)
    mu3_nom<-c(0,0,0)
    mu3_denom<-c(0,0,0)
    mu4_nom<-c(0,0,0)
    mu4_denom<-c(0,0,0)
    mu5_nom<-c(0,0,0)
    mu5_denom<-c(0,0,0)
    mu6_nom<-c(0,0,0)
    mu6_denom<-c(0,0,0)
    mu7_nom<-c(0,0,0)
    mu7_denom<-c(0,0,0)
    mu8_nom<-c(0,0,0)
    mu8_denom<-c(0,0,0)
    for (g in 1:length(g_clusternum)){
      mu1_nom<-mu1_nom+zi[g,1]*mg[g]*mug[g,]%*%solve(covh[,,1])
      mu1_denom<-mu1_denom+zi[g,1]*mg[g]*solve(covh[,,1])

      mu2_nom<-mu2_nom+zi[g,2]*mg[g]*mug[g,]%*%solve(covh[,,2])
      mu2_denom<-mu2_denom+zi[g,2]*mg[g]*solve(covh[,,2])

      mu3_nom<-mu3_nom+zi[g,3]*mg[g]*mug[g,]%*%solve(covh[,,3])
      mu3_denom<-mu3_denom+zi[g,3]*mg[g]*solve(covh[,,3])

      mu4_nom<-mu4_nom+zi[g,4]*mg[g]*mug[g,]%*%solve(covh[,,4])
      mu4_denom<-mu4_denom+zi[g,4]*mg[g]*solve(covh[,,4])

      mu5_nom<-mu5_nom+zi[g,5]*mg[g]*mug[g,]%*%solve(covh[,,5])
      mu5_denom<-mu5_denom+zi[g,5]*mg[g]*solve(covh[,,5])

      mu6_nom<-mu6_nom+zi[g,6]*mg[g]*mug[g,]%*%solve(covh[,,6])
      mu6_denom<-mu6_denom+zi[g,6]*mg[g]*solve(covh[,,6])

      mu7_nom<-mu7_nom+zi[g,7]*mg[g]*mug[g,]%*%solve(covh[,,7])
      mu7_denom<-mu7_denom+zi[g,7]*mg[g]*solve(covh[,,7])

      mu8_nom<-mu8_nom+zi[g,8]*mg[g]*mug[g,]%*%solve(covh[,,8])
      mu8_denom<-mu8_denom+zi[g,8]*mg[g]*solve(covh[,,8])
    }
    mu1_nom<-mu1_nom+2*lambdas[1]*neg_assum+2*lambdas[2]*(coefs[1]+coefs[2]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-muh[5,])+2*lambdas[3]*(coefs[1]+coefs[3]-1)*(coefs[1]*muh[2,]+coefs[3]*muh[4,]-muh[6,])+2*lambdas[4]*(coefs[2]+coefs[3]-1)*(coefs[2]*muh[3,]+coefs[3]*muh[4,]-muh[7,])+2*lambdas[5]*(coefs[1]+coefs[2]+coefs[3]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-muh[8,])
    mu1_denom<-mu1_denom+2*lambdas[1]*diag(1,nrow=dim_data)+2*lambdas[2]*(coefs[1]+coefs[2]-1)^2*diag(1,nrow=dim_data)+2*lambdas[3]*(coefs[1]+coefs[3]-1)^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[2]+coefs[3]-1)^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[1]+coefs[2]+coefs[3]-1)^2*diag(1,nrow=dim_data)

    mu2_nom<-mu2_nom+2*lambdas[2]*coefs[1]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[5,]-coefs[2]*muh[3,])+2*lambdas[3]*coefs[1]*((coefs[1]+coefs[3]-1)*muh[1,]+muh[6,]-coefs[3]*muh[4,])+2*lambdas[5]*coefs[1]*(muh[8,]-coefs[2]*muh[3,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu2_denom<-mu2_denom+2*lambdas[2]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[3]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[1])^2*diag(1,nrow=dim_data)

    mu3_nom<-mu3_nom+2*lambdas[2]*coefs[2]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[5,]-coefs[1]*muh[2,])+2*lambdas[4]*coefs[2]*((coefs[2]+coefs[3]-1)*muh[1,]+muh[7,]-coefs[3]*muh[4,])+2*lambdas[5]*coefs[2]*(muh[8,]-coefs[1]*muh[2,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu3_denom<-mu3_denom+2*lambdas[2]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[2])^2*diag(1,nrow=dim_data)

    mu4_nom<-mu4_nom+2*lambdas[3]*coefs[3]*((coefs[1]+coefs[3]-1)*muh[1,]+muh[6,]-coefs[1]*muh[2,])+2*lambdas[4]*coefs[3]*((coefs[2]+coefs[3]-1)*muh[1,]+muh[7,]-coefs[2]*muh[3,])+2*lambdas[5]*coefs[3]*(muh[8,]-coefs[1]*muh[2,]-coefs[2]*muh[3,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu4_denom<-mu4_denom+2*lambdas[3]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[3])^2*diag(1,nrow=dim_data)

    mu5_nom<-mu5_nom+2*lambdas[2]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-(coefs[1]+coefs[2]-1)*muh[1,])
    mu5_denom<-mu5_denom+2*lambdas[2]*diag(1,nrow=dim_data)

    mu6_nom<-mu6_nom+2*lambdas[3]*(coefs[1]*muh[2,]+coefs[3]*muh[4,]-(coefs[1]+coefs[3]-1)*muh[1,])
    mu6_denom<-mu6_denom+2*lambdas[3]*diag(1,nrow=dim_data)

    mu7_nom<-mu7_nom+2*lambdas[4]*(coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[2]+coefs[3]-1)*muh[1,])
    mu7_denom<-mu7_denom+2*lambdas[4]*diag(1,nrow=dim_data)

    mu8_nom<-mu8_nom+2*lambdas[5]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu8_denom<-mu8_denom+2*lambdas[5]*diag(1,nrow=dim_data)

    muh[1,]=mu1_nom%*%solve(mu1_denom)
    muh[2,]=mu2_nom%*%solve(mu2_denom)
    muh[3,]=mu3_nom%*%solve(mu3_denom)
    muh[4,]=mu4_nom%*%solve(mu4_denom)
    muh[5,]=mu5_nom%*%solve(mu5_denom)
    muh[6,]=mu6_nom%*%solve(mu6_denom)
    muh[7,]=mu7_nom%*%solve(mu7_denom)
    muh[8,]=mu8_nom%*%solve(mu8_denom)

    ## for Sigmas
    sigma_nom<-matrix(rep(0,dim_data^2),nrow=dim_data)
    sigma_denom<-0
    for (k in 1:cluster_num) {
      for (g in 1:length(g_clusternum)){
        sigma_nom<-sigma_nom+zi[g,k]*mg[g]*covg[,,g]+zi[g,k]*mg[g]*(mug[g,]-muh[k,])%*%(t(mug[g,]-muh[k,]))
        sigma_denom<-sigma_denom+zi[g,k]*mg[g]
      }
      covh[,, k] <- sigma_nom/sigma_denom
    }

    # change the likelihood
    change_ests<-c(change_ests,sum(apply(zi,1,function(x) log(sum(x)))))
    ## stopping criteria
    if(j>1){

      if(sum(abs(muh-old_mu))<eps){break}
    }

    old_mu <- muh


  }

  return(list(zi,muh,pih))
}

#' Internal Function 6
#'
#' This function merges 4-d dPCR data
#' @param data A matrix or data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @param eps the convergence threshold
#' @param max_iter maximum number of iterations
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A list of membership probability, cluster center, merging probability
#' @keywords internal
HMM_merge_4d<-function(data,cluster_num,fp,eps=10^(-10),max_iter=1000,lambdas=rep(2,12),coefs=rep(1,4)) {
  change_ests=NULL
  data<-as.matrix(data)
  dim_data<-ncol(data)

  pih = rep(1/cluster_num,cluster_num)
  g_clusternum<-unique(fp$peaks.cluster)
  mg<-table(fp$peaks.cluster)
  mug<-fp$peaks$mu
  covg <- array(0,dim=c(dim_data,dim_data,length(g_clusternum)))
  for (i in 1:length(g_clusternum)){
    # tryCatch( {
    #   covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    #   print(covg[,,g_clusternum[i]])
    # }, error = function(e) {
    #   covg[,,g_clusternum[i]] <<- covg[,,g_clusternum[i-1]]
    #   print(covg[,,g_clusternum[i-1]])
    # })
    tryCatch( {
      covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    }, error = function(e) {
      covg[,,g_clusternum[i]] <<- diag(0.001,nrow=dim_data)
    })
  }

  #initialization
  ## find the negative population
  min_val<-apply(data,2,min)
  dist_orig<-apply(mug,1,function(x) {sqrt(sum((x-min_val)^2))})
  neg_assum<-mug[which.min(dist_orig),]

  muh<-matrix(0,nrow=cluster_num,ncol=dim_data)
  muh[1,]<-neg_assum
  muh[2,]<-c(0.75*(min(mug[,1])+max(mug[,1])),neg_assum[2],neg_assum[3],neg_assum[4])
  muh[3,]<-c(neg_assum[1],0.75*(min(mug[,2])+max(mug[,2])),neg_assum[3],neg_assum[4])
  muh[4,]<-c(neg_assum[1],neg_assum[2],0.75*(min(mug[,3])+max(mug[,3])),neg_assum[4])
  muh[5,]<-c(neg_assum[1],neg_assum[2],neg_assum[3],0.75*(min(mug[,4])+max(mug[,4])))

  muh[6,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]-(coefs[1]+coefs[2]-1)*muh[1,]
  muh[7,]<-coefs[1]*muh[2,]+coefs[3]*muh[4,]-(coefs[1]+coefs[3]-1)*muh[1,]
  muh[8,]<-coefs[1]*muh[2,]+coefs[4]*muh[5,]-(coefs[1]+coefs[4]-1)*muh[1,]
  muh[9,]<-coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[2]+coefs[3]-1)*muh[1,]
  muh[10,]<-coefs[2]*muh[3,]+coefs[4]*muh[5,]-(coefs[2]+coefs[4]-1)*muh[1,]
  muh[11,]<-coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[3]+coefs[4]-1)*muh[1,]

  muh[12,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,]
  muh[13,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[4]*muh[5,]-(coefs[1]+coefs[2]+coefs[4]-1)*muh[1,]
  muh[14,]<-coefs[1]*muh[2,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[1]+coefs[3]+coefs[4]-1)*muh[1,]
  muh[15,]<-coefs[2]*muh[3,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[2]+coefs[3]+coefs[4]-1)*muh[1,]

  muh[16,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*muh[1,]


  covh <- array(cov(mug), dim = c(dim_data, dim_data, cluster_num))

  zi<-matrix(0,nrow=length(g_clusternum), ncol=cluster_num)

  ## start the EM algorithm

  for (j in 1:max_iter){
    if (j>=max_iter){
      print("Warning: the algorithm fails to converge")
    }
    ## E step:
    for (g in 1:length(g_clusternum)){
      for (k in 1:cluster_num){
        inv_zi_tmp<-0
        for (l in 1:cluster_num){
          tmp=exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
          inv_zi_tmp=inv_zi_tmp+exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
        }
        zi[g,k]=1/inv_zi_tmp
      }
    }


    ## M step:
    pih = apply(zi,2,sum)/length(g_clusternum)
    pih[which(pih==0)]<-(10^(-10))
    ## for the mus
    mus_nom<-matrix(0,nrow=cluster_num,ncol=dim_data)
    mus_denom<-array(0, dim = c(dim_data, dim_data, cluster_num))


    for (g in 1:length(g_clusternum)){
      for (h in 1:cluster_num){
        mus_nom[h,]<-mus_nom[h,]+zi[g,h]*mg[g]*mug[g,]%*%solve(covh[,,h])
        mus_denom[,,h]<-mus_denom[,,h]+zi[g,h]*mg[g]*solve(covh[,,h])
      }
    }
    mus_nom[1,]<-mus_nom[1,]+2*lambdas[1]*neg_assum+2*lambdas[2]*(coefs[1]+coefs[2]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-muh[6,])+2*lambdas[3]*(coefs[1]+coefs[3]-1)*(coefs[1]*muh[2,]+coefs[3]*muh[4,]-muh[7,])+2*lambdas[4]*(coefs[1]+coefs[4]-1)*(coefs[1]*muh[2,]+coefs[4]*muh[5,]-muh[8,])
    +2*lambdas[5]*(coefs[2]+coefs[3]-1)*(coefs[2]*muh[3,]+coefs[3]*muh[4,]-muh[9,])+2*lambdas[6]*(coefs[2]+coefs[4]-1)*(coefs[2]*muh[3,]+coefs[4]*muh[5,]-muh[10,])+2*lambdas[7]*(coefs[3]+coefs[4]-1)*(coefs[3]*muh[4,]+coefs[4]*muh[5,]-muh[11,])+2*lambdas[8]*(coefs[1]+coefs[2]+coefs[3]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-muh[12,])
    +2*lambdas[9]*(coefs[1]+coefs[2]+coefs[4]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[4]*muh[5,]-muh[13,])+2*lambdas[10]*(coefs[1]+coefs[3]+coefs[4]-1)*(coefs[1]*muh[2,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-muh[14,])+2*lambdas[11]*(coefs[2]+coefs[3]+coefs[4]-1)*(coefs[2]*muh[3,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-muh[15,])
    +2*lambdas[12]*(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-muh[16,])
    mus_denom[,,1]<-mus_denom[,,1]+2*lambdas[1]*diag(1,nrow=dim_data)+2*lambdas[2]*(coefs[1]+coefs[2]-1)^2*diag(1,nrow=dim_data)+2*lambdas[3]*(coefs[1]+coefs[3]-1)^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[1]+coefs[4]-1)^2*diag(1,nrow=dim_data)
    +2*lambdas[5]*(coefs[2]+coefs[3]-1)^2*diag(1,nrow=dim_data)+2*lambdas[6]*(coefs[2]+coefs[4]-1)^2*diag(1,nrow=dim_data)+2*lambdas[7]*(coefs[3]+coefs[4]-1)^2*diag(1,nrow=dim_data)+2*lambdas[8]*(coefs[1]+coefs[2]+coefs[3]-1)^2*diag(1,nrow=dim_data)
    +2*lambdas[9]*(coefs[1]+coefs[2]+coefs[4]-1)^2*diag(1,nrow=dim_data)+2*lambdas[10]*(coefs[1]+coefs[3]+coefs[4]-1)^2*diag(1,nrow=dim_data)+2*lambdas[11]*(coefs[2]+coefs[3]+coefs[4]-1)^2*diag(1,nrow=dim_data)+2*lambdas[12]*(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)^2*diag(1,nrow=dim_data)

    mus_nom[2,]<-mus_nom[2,]+2*lambdas[2]*coefs[1]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[6,]-coefs[2]*muh[3,])+2*lambdas[3]*coefs[1]*((coefs[1]+coefs[3]-1)*muh[1,]+muh[7,]-coefs[3]*muh[4,])+2*lambdas[4]*coefs[1]*((coefs[1]+coefs[4]-1)*muh[1,]+muh[8,]-coefs[4]*muh[5,])
    +2*lambdas[8]*coefs[1]*(muh[12,]-coefs[2]*muh[3,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])+2*lambdas[9]*coefs[1]*(muh[13,]-coefs[2]*muh[3,]-coefs[4]*muh[5,]+(coefs[1]+coefs[2]+coefs[4]-1)*muh[1,])+2*lambdas[10]*coefs[1]*(muh[14,]-coefs[3]*muh[4,]-coefs[4]*muh[5,]+(coefs[1]+coefs[3]+coefs[4]-1)*muh[1,])
    +2*lambdas[12]*coefs[1]*(muh[16,]-coefs[2]*muh[3,]-coefs[3]*muh[4,]-coefs[4]*muh[5,]+(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,2]<- mus_denom[,,2]+2*lambdas[2]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[3]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[1])^2*diag(1,nrow=dim_data)
    +2*lambdas[8]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[9]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[10]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[12]*(coefs[1])^2*diag(1,nrow=dim_data)

    mus_nom[3,]<-mus_nom[3,]+2*lambdas[2]*coefs[2]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[6,]-coefs[1]*muh[2,])+2*lambdas[5]*coefs[2]*((coefs[2]+coefs[3]-1)*muh[1,]+muh[9,]-coefs[3]*muh[4,])+2*lambdas[6]*coefs[2]*((coefs[2]+coefs[4]-1)*muh[1,]+muh[10,]-coefs[4]*muh[5,])
    +2*lambdas[8]*coefs[2]*(muh[12,]-coefs[1]*muh[2,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])+2*lambdas[9]*coefs[2]*(muh[13,]-coefs[1]*muh[2,]-coefs[4]*muh[5,]+(coefs[1]+coefs[2]+coefs[4]-1)*muh[1,])+2*lambdas[11]*coefs[2]*(muh[15,]-coefs[3]*muh[4,]-coefs[4]*muh[5,]+(coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    +2*lambdas[12]*coefs[2]*(muh[16,]-coefs[1]*muh[2,]-coefs[3]*muh[4,]-coefs[4]*muh[5,]+(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,3]<- mus_denom[,,3]+2*lambdas[2]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[6]*(coefs[2])^2*diag(1,nrow=dim_data)
    +2*lambdas[8]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[9]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[11]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[12]*(coefs[2])^2*diag(1,nrow=dim_data)

    mus_nom[4,]<-mus_nom[4,]+2*lambdas[3]*coefs[3]*((coefs[1]+coefs[3]-1)*muh[1,]+muh[7,]-coefs[1]*muh[2,])+2*lambdas[5]*coefs[3]*((coefs[2]+coefs[3]-1)*muh[1,]+muh[9,]-coefs[2]*muh[3,])+2*lambdas[7]*coefs[3]*((coefs[3]+coefs[4]-1)*muh[1,]+muh[11,]-coefs[4]*muh[5,])
    +2*lambdas[8]*coefs[3]*(muh[12,]-coefs[1]*muh[2,]-coefs[2]*muh[3,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])+2*lambdas[10]*coefs[3]*(muh[14,]-coefs[1]*muh[2,]-coefs[4]*muh[5,]+(coefs[1]+coefs[3]+coefs[4]-1)*muh[1,])+2*lambdas[11]*coefs[3]*(muh[15,]-coefs[2]*muh[3,]-coefs[4]*muh[5,]+(coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    +2*lambdas[12]*coefs[3]*(muh[16,]-coefs[1]*muh[2,]-coefs[2]*muh[3,]-coefs[4]*muh[5,]+(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,4]<- mus_denom[,,4]+2*lambdas[3]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[7]*(coefs[3])^2*diag(1,nrow=dim_data)
    +2*lambdas[8]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[10]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[11]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[12]*(coefs[3])^2*diag(1,nrow=dim_data)

    mus_nom[5,]<-mus_nom[5,]+2*lambdas[4]*coefs[4]*((coefs[1]+coefs[4]-1)*muh[1,]+muh[8,]-coefs[1]*muh[2,])+2*lambdas[6]*coefs[4]*((coefs[2]+coefs[4]-1)*muh[1,]+muh[10,]-coefs[2]*muh[3,])+2*lambdas[7]*coefs[4]*((coefs[3]+coefs[4]-1)*muh[1,]+muh[11,]-coefs[3]*muh[4,])
    +2*lambdas[9]*coefs[4]*(muh[13,]-coefs[1]*muh[2,]-coefs[2]*muh[3,]+(coefs[1]+coefs[2]+coefs[4]-1)*muh[1,])+2*lambdas[10]*coefs[4]*(muh[14,]-coefs[1]*muh[2,]-coefs[3]*muh[4,]+(coefs[1]+coefs[3]+coefs[4]-1)*muh[1,])+2*lambdas[11]*coefs[4]*(muh[15,]-coefs[2]*muh[3,]-coefs[3]*muh[4,]+(coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    +2*lambdas[12]*coefs[4]*(muh[16,]-coefs[1]*muh[2,]-coefs[2]*muh[3,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,5]<- mus_denom[,,5]+2*lambdas[4]*(coefs[4])^2*diag(1,nrow=dim_data)+2*lambdas[6]*(coefs[4])^2*diag(1,nrow=dim_data)+2*lambdas[7]*(coefs[4])^2*diag(1,nrow=dim_data)
    +2*lambdas[9]*(coefs[4])^2*diag(1,nrow=dim_data)+2*lambdas[10]*(coefs[4])^2*diag(1,nrow=dim_data)+2*lambdas[11]*(coefs[4])^2*diag(1,nrow=dim_data)+2*lambdas[12]*(coefs[4])^2*diag(1,nrow=dim_data)

    mus_nom[6,]<-mus_nom[6,]+2*lambdas[2]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-(coefs[1]+coefs[2]-1)*muh[1,])
    mus_denom[,,6]<- mus_denom[,,6]+2*lambdas[2]*diag(1,nrow=dim_data)

    mus_nom[7,]<-mus_nom[7,]+2*lambdas[3]*(coefs[1]*muh[2,]+coefs[3]*muh[4,]-(coefs[1]+coefs[3]-1)*muh[1,])
    mus_denom[,,7]<- mus_denom[,,7]+2*lambdas[3]*diag(1,nrow=dim_data)

    mus_nom[8,]<-mus_nom[8,]+2*lambdas[4]*(coefs[1]*muh[2,]+coefs[4]*muh[5,]-(coefs[1]+coefs[4]-1)*muh[1,])
    mus_denom[,,8]<- mus_denom[,,8]+2*lambdas[4]*diag(1,nrow=dim_data)

    mus_nom[9,]<-mus_nom[9,]+2*lambdas[5]*(coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[2]+coefs[3]-1)*muh[1,])
    mus_denom[,,9]<- mus_denom[,,9]+2*lambdas[5]*diag(1,nrow=dim_data)

    mus_nom[10,]<-mus_nom[10,]+2*lambdas[6]*(coefs[2]*muh[3,]+coefs[4]*muh[5,]-(coefs[2]+coefs[4]-1)*muh[1,])
    mus_denom[,,10]<- mus_denom[,,10]+2*lambdas[6]*diag(1,nrow=dim_data)

    mus_nom[11,]<-mus_nom[11,]+2*lambdas[7]*(coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,11]<- mus_denom[,,11]+2*lambdas[7]*diag(1,nrow=dim_data)

    mus_nom[12,]<-mus_nom[12,]+2*lambdas[8]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mus_denom[,,12]<- mus_denom[,,12]+2*lambdas[8]*diag(1,nrow=dim_data)

    mus_nom[13,]<-mus_nom[13,]+2*lambdas[9]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[4]*muh[5,]-(coefs[1]+coefs[2]+coefs[4]-1)*muh[1,])
    mus_denom[,,13]<- mus_denom[,,13]+2*lambdas[9]*diag(1,nrow=dim_data)

    mus_nom[14,]<-mus_nom[14,]+2*lambdas[10]*(coefs[1]*muh[2,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[1]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,14]<- mus_denom[,,14]+2*lambdas[10]*diag(1,nrow=dim_data)

    mus_nom[15,]<-mus_nom[15,]+2*lambdas[11]*(coefs[2]*muh[3,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,15]<- mus_denom[,,15]+2*lambdas[11]*diag(1,nrow=dim_data)

    mus_nom[16,]<-mus_nom[16,]+2*lambdas[12]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]+coefs[4]*muh[5,]-(coefs[1]+coefs[2]+coefs[3]+coefs[4]-1)*muh[1,])
    mus_denom[,,16]<- mus_denom[,,16]+2*lambdas[12]*diag(1,nrow=dim_data)

    for (h in 1:cluster_num){
      muh[h,]=mus_nom[h,]%*%solve(mus_denom[,,h])
    }



    ## for Sigmas
    sigma_nom<-matrix(rep(0,dim_data^2),nrow=dim_data)
    sigma_denom<-0
    for (k in 1:cluster_num) {
      for (g in 1:length(g_clusternum)){
        sigma_nom<-sigma_nom+zi[g,k]*mg[g]*covg[,,g]+zi[g,k]*mg[g]*(mug[g,]-muh[k,])%*%(t(mug[g,]-muh[k,]))
        sigma_denom<-sigma_denom+zi[g,k]*mg[g]
      }
      covh[,, k] <- sigma_nom/sigma_denom
    }

    # change the likelihood
    change_ests<-c(change_ests,sum(apply(zi,1,function(x) log(sum(x)))))
    ## stopping criteria
    if(j>1){

      if(sum(abs(muh-old_mu))<eps){break}
    }

    old_mu <- muh


  }

  return(list(zi,muh,pih))
}

#' Internal Function 7
#'
#' This function merges 2-d and 3-target dPCR data
#' @param data A matrix or data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @param eps the convergence threshold
#' @param max_iter maximum number of iterations
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A list of membership probability, cluster center, merging probability
#' @keywords internal
HMM_merge_higher_order<-function(data,cluster_num,fp,eps=10^(-10),max_iter=1000,lambdas=rep(2,5),coefs=rep(1,3)) {
  change_ests=NULL
  data<-as.matrix(data)
  dim_data<-ncol(data)

  pih = rep(1/cluster_num,cluster_num)
  g_clusternum<-unique(fp$peaks.cluster)
  mg<-table(fp$peaks.cluster)
  mug<-fp$peaks$mu
  covg <- array(0,dim=c(dim_data,dim_data,length(g_clusternum)))
  for (i in 1:length(g_clusternum)){
    # tryCatch( {
    #   covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    #   print(covg[,,g_clusternum[i]])
    # }, error = function(e) {
    #   covg[,,g_clusternum[i]] <<- covg[,,g_clusternum[i-1]]
    #   print(covg[,,g_clusternum[i-1]])
    # })
    tryCatch( {
      covg[,,g_clusternum[i]] <- cov(data[fp$peaks.cluster==g_clusternum[i],])
    }, error = function(e) {
      covg[,,g_clusternum[i]] <<- diag(0.001,nrow=dim_data)
    })
  }

  #initialization
  ## find the negative population
  min_val<-apply(data,2,min)
  dist_orig<-apply(mug,1,function(x) {sqrt(sum((x-min_val)^2))})
  neg_assum<-mug[which.min(dist_orig),]

  muh<-matrix(0,nrow=cluster_num,ncol=dim_data)
  muh[1,]<-neg_assum
  muh[2,]<-c(0.25*min(mug[,1])+0.75*max(mug[,1]),neg_assum[2])
  muh[3,]<-c(neg_assum[1],0.25*min(mug[,2])+0.75*max(mug[,2]))
  muh[4,]<-c(0.5*(neg_assum[1]+muh[2,1]),0.5*(neg_assum[2]+muh[3,2]))
  muh[5,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+(1-coefs[1]-coefs[2])*muh[1,]
  muh[6,]<-coefs[1]*muh[2,]+coefs[3]*muh[4,]+(1-coefs[1]-coefs[3])*muh[1,]
  muh[7,]<-coefs[2]*muh[3,]+coefs[3]*muh[4,]+(1-coefs[2]-coefs[3])*muh[1,]
  muh[8,]<-coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]+(1-coefs[1]-coefs[2]-coefs[3])*muh[1,]

  covh <- array(cov(mug), dim = c(dim_data, dim_data, cluster_num))

  zi<-matrix(0,nrow=length(g_clusternum), ncol=cluster_num)

  ## start the EM algorithm

  for (j in 1:max_iter){
    if (j>=max_iter){
      print("Warning: the algorithm fails to converge")
    }
    ## E step:
    for (g in 1:length(g_clusternum)){
      for (k in 1:cluster_num){
        inv_zi_tmp<-0
        for (l in 1:cluster_num){
          tmp=exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
          inv_zi_tmp=inv_zi_tmp+exp(log(pih[l])-log(pih[k])+mg[g]*(dmvnorm(t(mug[g,]),t(muh[l,]),covh[,,l],log=TRUE)-0.5*tr(solve(covh[,,l])%*%covg[,,g])-dmvnorm(t(mug[g,]),t(muh[k,]),covh[,,k],log=TRUE)+0.5*tr(solve(covh[,,k])%*%covg[,,g])))
        }
        zi[g,k]=1/inv_zi_tmp
      }
    }


    ## M step:
    pih = apply(zi,2,sum)/length(g_clusternum)
    pih[which(pih==0)]<-(10^(-10))
    ## for the mus
    mu1_nom<-c(0,0)
    mu1_denom<-c(0,0)
    mu2_nom<-c(0,0)
    mu2_denom<-c(0,0)
    mu3_nom<-c(0,0)
    mu3_denom<-c(0,0)
    mu4_nom<-c(0,0)
    mu4_denom<-c(0,0)
    mu5_nom<-c(0,0)
    mu5_denom<-c(0,0)
    mu6_nom<-c(0,0)
    mu6_denom<-c(0,0)
    mu7_nom<-c(0,0)
    mu7_denom<-c(0,0)
    mu8_nom<-c(0,0)
    mu8_denom<-c(0,0)
    for (g in 1:length(g_clusternum)){
      mu1_nom<-mu1_nom+zi[g,1]*mg[g]*mug[g,]%*%solve(covh[,,1])
      mu1_denom<-mu1_denom+zi[g,1]*mg[g]*solve(covh[,,1])

      mu2_nom<-mu2_nom+zi[g,2]*mg[g]*mug[g,]%*%solve(covh[,,2])
      mu2_denom<-mu2_denom+zi[g,2]*mg[g]*solve(covh[,,2])

      mu3_nom<-mu3_nom+zi[g,3]*mg[g]*mug[g,]%*%solve(covh[,,3])
      mu3_denom<-mu3_denom+zi[g,3]*mg[g]*solve(covh[,,3])

      mu4_nom<-mu4_nom+zi[g,4]*mg[g]*mug[g,]%*%solve(covh[,,4])
      mu4_denom<-mu4_denom+zi[g,4]*mg[g]*solve(covh[,,4])

      mu5_nom<-mu5_nom+zi[g,5]*mg[g]*mug[g,]%*%solve(covh[,,5])
      mu5_denom<-mu5_denom+zi[g,5]*mg[g]*solve(covh[,,5])

      mu6_nom<-mu6_nom+zi[g,6]*mg[g]*mug[g,]%*%solve(covh[,,6])
      mu6_denom<-mu6_denom+zi[g,6]*mg[g]*solve(covh[,,6])

      mu7_nom<-mu7_nom+zi[g,7]*mg[g]*mug[g,]%*%solve(covh[,,7])
      mu7_denom<-mu7_denom+zi[g,7]*mg[g]*solve(covh[,,7])

      mu8_nom<-mu8_nom+zi[g,8]*mg[g]*mug[g,]%*%solve(covh[,,8])
      mu8_denom<-mu8_denom+zi[g,8]*mg[g]*solve(covh[,,8])
    }
    mu1_nom<-mu1_nom+2*lambdas[1]*neg_assum+2*lambdas[2]*(coefs[1]+coefs[2]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-muh[5,])+2*lambdas[3]*(coefs[1]+coefs[3]-1)*(coefs[1]*muh[2,]+coefs[3]*muh[4,]-muh[6,])+2*lambdas[4]*(coefs[2]+coefs[3]-1)*(coefs[2]*muh[3,]+coefs[3]*muh[4,]-muh[7,])+2*lambdas[5]*(coefs[1]+coefs[2]+coefs[3]-1)*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-muh[8,])
    mu1_denom<-mu1_denom+2*lambdas[1]*diag(1,nrow=dim_data)+2*lambdas[2]*(coefs[1]+coefs[2]-1)^2*diag(1,nrow=dim_data)+2*lambdas[3]*(coefs[1]+coefs[3]-1)^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[2]+coefs[3]-1)^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[1]+coefs[2]+coefs[3]-1)^2*diag(1,nrow=dim_data)

    mu2_nom<-mu2_nom+2*lambdas[2]*coefs[1]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[5,]-coefs[2]*muh[3,])+2*lambdas[3]*coefs[1]*((coefs[1]+coefs[3]-1)*muh[1,]+muh[6,]-coefs[3]*muh[4,])+2*lambdas[5]*coefs[1]*(muh[8,]-coefs[2]*muh[3,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu2_denom<-mu2_denom+2*lambdas[2]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[3]*(coefs[1])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[1])^2*diag(1,nrow=dim_data)

    mu3_nom<-mu3_nom+2*lambdas[2]*coefs[2]*((coefs[1]+coefs[2]-1)*muh[1,]+muh[5,]-coefs[1]*muh[2,])+2*lambdas[4]*coefs[2]*((coefs[2]+coefs[3]-1)*muh[1,]+muh[7,]-coefs[3]*muh[4,])+2*lambdas[5]*coefs[2]*(muh[8,]-coefs[1]*muh[2,]-coefs[3]*muh[4,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu3_denom<-mu3_denom+2*lambdas[2]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[2])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[2])^2*diag(1,nrow=dim_data)

    mu4_nom<-mu4_nom+2*lambdas[3]*coefs[3]*((coefs[1]+coefs[3]-1)*muh[1,]+muh[6,]-coefs[1]*muh[2,])+2*lambdas[4]*coefs[3]*((coefs[2]+coefs[3]-1)*muh[1,]+muh[7,]-coefs[2]*muh[3,])+2*lambdas[5]*coefs[3]*(muh[8,]-coefs[1]*muh[2,]-coefs[2]*muh[3,]+(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu4_denom<-mu4_denom+2*lambdas[3]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[4]*(coefs[3])^2*diag(1,nrow=dim_data)+2*lambdas[5]*(coefs[3])^2*diag(1,nrow=dim_data)

    mu5_nom<-mu5_nom+2*lambdas[2]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]-(coefs[1]+coefs[2]-1)*muh[1,])
    mu5_denom<-mu5_denom+2*lambdas[2]*diag(1,nrow=dim_data)

    mu6_nom<-mu6_nom+2*lambdas[3]*(coefs[1]*muh[2,]+coefs[3]*muh[4,]-(coefs[1]+coefs[3]-1)*muh[1,])
    mu6_denom<-mu6_denom+2*lambdas[3]*diag(1,nrow=dim_data)

    mu7_nom<-mu7_nom+2*lambdas[4]*(coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[2]+coefs[3]-1)*muh[1,])
    mu7_denom<-mu7_denom+2*lambdas[4]*diag(1,nrow=dim_data)

    mu8_nom<-mu8_nom+2*lambdas[5]*(coefs[1]*muh[2,]+coefs[2]*muh[3,]+coefs[3]*muh[4,]-(coefs[1]+coefs[2]+coefs[3]-1)*muh[1,])
    mu8_denom<-mu8_denom+2*lambdas[5]*diag(1,nrow=dim_data)

    muh[1,]=mu1_nom%*%solve(mu1_denom)
    muh[2,]=mu2_nom%*%solve(mu2_denom)
    muh[3,]=mu3_nom%*%solve(mu3_denom)
    muh[4,]=mu4_nom%*%solve(mu4_denom)
    muh[5,]=mu5_nom%*%solve(mu5_denom)
    muh[6,]=mu6_nom%*%solve(mu6_denom)
    muh[7,]=mu7_nom%*%solve(mu7_denom)
    muh[8,]=mu8_nom%*%solve(mu8_denom)

    ## for Sigmas
    sigma_nom<-matrix(rep(0,dim_data^2),nrow=dim_data)
    sigma_denom<-0
    for (k in 1:cluster_num) {
      for (g in 1:length(g_clusternum)){
        sigma_nom<-sigma_nom+zi[g,k]*mg[g]*covg[,,g]+zi[g,k]*mg[g]*(mug[g,]-muh[k,])%*%(t(mug[g,]-muh[k,]))
        sigma_denom<-sigma_denom+zi[g,k]*mg[g]
      }
      covh[,, k] <- sigma_nom/sigma_denom
    }

    # change the likelihood
    change_ests<-c(change_ests,sum(apply(zi,1,function(x) log(sum(x)))))
    ## stopping criteria
    if(j>1){

      if(sum(abs(muh-old_mu))<eps){break}
    }

    old_mu <- muh


  }

  return(list(zi,muh,pih))
}
