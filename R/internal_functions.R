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
    
    clust_pos<- data %>% group_by(cluster) %>% summarise(mean_x=mean(channel1),mean_y=mean(channel2))
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
fp_search <- function(data, cluster_num = 16) {
    # Simple 2D objective function
    obj.fun.sil <- makeSingleObjectiveFunction(
        fn = function(pars) {
            fp <- flowPeaks(data, tol = pars[1], h0 = pars[2], h = pars[3])
            deviation <- cluster_num + 2 - length(unique(fp$peaks.cluster))
            fp_sil_optimval <- tryCatch({
                deviation^2
            }, error = function(e) {
                10
            })
            return(fp_sil_optimval)
        },
        par.set = makeParamSet(
            makeNumericParam("x1", lower = 0, upper = 1),
            makeNumericParam("x2", lower = 0.1, upper = 5),
            makeNumericParam("x3", lower = 0.1, upper = 5)
        )
    )
    
    # Create base control object
    ctrl <- makeMBOControl()
    # Set MBO control termination
    ctrl <- setMBOControlTermination(ctrl, iters = 4L)
    # Use 500 points in the focus search (should be sufficient for 2D)
    ctrl <- setMBOControlInfill(ctrl, opt.focussearch.points = 1000)
    # Create initial design
    des <- generateDesign(n = 15L, getParamSet(obj.fun.sil), fun = lhs::maximinLHS)
    # Start MBO
    res <- mbo(obj.fun.sil, design = des, control = ctrl)
    bo_sil_result <- unlist(res$x)
    return(bo_sil_result)
}

#' Internal Function 4
#'
#' This function outputs vectors and weights that will be used in EM algorithm
#' @param coefs coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @param mus The cluster centers of primary targets
#' @param cluster_num The expected maximum number of clusters.
#' @param dim_data dimension of the dataset
#' @return A list of vectors and weights
#' @keywords internal
combined_vectors<-function(coefs,mus,cluster_num,dim_data){
    comb_num<-cluster_num-log2(cluster_num)-1
    primary_tar<-log2(cluster_num)
    comb_seq<-seq(2,primary_tar)
    
    coef_matrix<-matrix(0,nrow=primary_tar+1,ncol=primary_tar+1)
    coef_matrix[,1]<-c(1,coefs)
    diag(coef_matrix)<-c(1,-coefs)
    all_combinations <- list()
    mat_coef<-NULL
    
    coef_tmp<-c(1,rep(0,primary_tar))
    for (i in comb_seq){
        comb_tmp<-combn(seq_len(primary_tar), i)
        all_combinations[[i]] <- comb_tmp
        comb_ncol<-ncol(comb_tmp)
        for (j in seq_len(comb_ncol)){
            mat_coef_tmp<-coef_tmp
            mat_coef_tmp[comb_tmp[,j]+1]<-(-1)
            mat_coef<-rbind(mat_coef,mat_coef_tmp)
        }
    }
    
    weights <- mat_coef%*%coef_matrix
    combined_vecs<-weights%*%mus
    
    return(list(combined_vecs,weights))
}



#' Internal Function 5
#'
#' This function intialize the parameters for the main clustering function
#' @param data A matrix or data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @param base_clust The results of base clustering
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A list of initial parameters for the EM algorithm
#' @keywords internal
GMM_init <- function(data, cluster_num, base_clust, coefs) {
    change_ests <- NULL
    data <- as.matrix(data)
    dim_data <- ncol(data)
    
    pih <- rep(1 / cluster_num, cluster_num)
    g_clusternum <- unique(base_clust$cluster)
    mg <- table(base_clust$cluster)
    mug <- base_clust$mu
    covg <- array(0, dim = c(dim_data, dim_data, length(g_clusternum)))
    
    for (i in seq_len(length(g_clusternum))) {
        covg[, , g_clusternum[i]] <- tryCatch({
            cov(data[base_clust$cluster == g_clusternum[i], ])
        }, error = function(e) {
            diag(0.001, nrow = dim_data)
        })
    }
    # Initialization
    min_val <- apply(data, 2, min)
    dist_orig <- apply(mug, 1, function(x) sqrt(sum((x - min_val)^2)))
    neg_assum <- mug[which.min(dist_orig), ]
    
    muh <- matrix(neg_assum, nrow = cluster_num, ncol = dim_data, byrow = TRUE)
    muh[2, 1] <- 0.75 * (min(mug[, 1]) + max(mug[, 1]))
    muh[3, 2] <-  0.75 * (min(mug[, 2]) + max(mug[, 2]))
    if (cluster_num == 8){
        if (dim_data==2){
            muh[4, ] <- c(0.5*(neg_assum[1]+muh[2,1]),0.5*(neg_assum[2]+muh[3,2]))
        } else if (dim_data==3) {
            muh[4, 3] <-  0.75 * (min(mug[, 3]) + max(mug[, 3]))
        }
    } else if (cluster_num >= 16) {
        muh[4, 3] <-  0.75 * (min(mug[, 3]) + max(mug[, 3]))
        muh[5, 4] <-  0.75 * (min(mug[, 4]) + max(mug[, 4]))
        if (cluster_num>=32) {
            muh[6, 5] <-  0.75 * (min(mug[, 5]) + max(mug[, 5]))
        }
        if (cluster_num == 64) {
            muh[7, 6] <-  0.75 * (min(mug[, 6]) + max(mug[, 6]))
        }
    }
    combined_results <- combined_vectors(coefs=coefs,mus=muh[seq_len(log2(cluster_num)+1),],cluster_num,dim_data)
    muh[-seq_len(log2(cluster_num)+1),] <- combined_results[[1]]
    weights <-combined_results[[2]]
    covh <- array(cov(mug), dim = c(dim_data, dim_data, cluster_num))
    
    return(list(pih,muh,covh,g_clusternum,mg,mug,covg, weights,neg_assum))
}


#' Internal Function 6
#'
#' This function compute the necessary elements for estep function 
#' @param g cluster index
#' @param k cluster index
#' @param cluster_num The expected maximum number of clusters
#' @param mg cluster sizes of base clustering result
#' @param log_pih log pih (the probability of cluster g belonging at level l+1 to cluster h at level l)
#' @param mug_t the transposed matrix of cluster centers at level l+1
#' @param muh_t the transposed matrix of cluster centers at level l
#' @param covg  the covariance matrix of clusters at level l+1
#' @param covh  the covariance matrix of clusters at level l
#' @return A vector of intermediate values for zi calculation in estep function
#' @keywords internal
compute_tmp_matrix <- function(g, k, cluster_num, mg, log_pih, mug_t, muh_t, covh, covg) {
    exp_diff <- numeric(cluster_num)
    
    for (l in seq_len(cluster_num)) {
        log_dmvnorm_l <- dmvnorm(mug_t[, g], muh_t[, l], covh[, , l], log = TRUE)
        log_dmvnorm_k <- dmvnorm(mug_t[, g], muh_t[, k], covh[, , k], log = TRUE)
        term_l <- -0.5 * tr(solve(covh[, , l]) %*% covg[, , g])
        term_k <- 0.5 * tr(solve(covh[, , k]) %*% covg[, , g])
        
        exp_diff[l] <- exp(log_pih[l] - log_pih[k] + mg[g] * (log_dmvnorm_l + term_l - log_dmvnorm_k + term_k))
    }
    
    return(exp_diff)
}

#' Internal Function 7
#'
#' This function calculates zi in E-step of EM algorithm
#' @param g_clusternum cluster labels from base clustering 
#' @param cluster_num The expected maximum number of clusters
#' @param pih the probability of cluster g belonging at level l+1 to cluster h at level l
#' @param muh the matrix of cluster centers at level l
#' @param covh the covariance matrix of clusters at level l
#' @param mg cluster sizes of base clustering result
#' @param mug the matrix of cluster centers at level l+1
#' @param covg the covariance matrix of clusters at level l+1
#' @return zi for estep in EM algorithm
#' @keywords internal
estep <- function(g_clusternum,cluster_num,pih,muh,covh,mg,mug,covg){
    log_pih <- log(pih)
    mug_t <- t(mug)
    muh_t <- t(muh)
    
    zi <- matrix(0, nrow = length(g_clusternum), ncol = cluster_num)
    # E step
    for (g in seq_len(length(g_clusternum))) {
        for (k in seq_len(cluster_num)) {
            tmp_matrix <- compute_tmp_matrix(g, k, cluster_num, mg, log_pih, mug_t, muh_t, covh, covg)
            inv_zi_tmp <- sum(tmp_matrix)
            zi[g, k] <- 1 / inv_zi_tmp
        }
    }
    return(zi)
}


#' Internal Function 8
#'
#' This function calculates mu in M-step of EM algorithm
#' @param zi the expected log-likelihood found on the E step
#' @param g_clusternum cluster labels from base clustering 
#' @param cluster_num The expected maximum number of clusters
#' @param dim_data the dimension of the dataset
#' @param weights combinations of coefficients of the cluster centers
#' @param muh the matrix of cluster centers at level l
#' @param covh the covariance matrix of clusters at level l
#' @param mg cluster sizes of base clustering result
#' @param mug the matrix of cluster centers at level l+1
#' @param neg_assum the estimated cluster center of negative population
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return muh the cluster centers at level l in the EM algorithm
#' @keywords internal
mstep_mu<-function(zi,g_clusternum,dim_data,cluster_num,weights,muh,covh,mg,mug,neg_assum,lambdas,coefs){
    # M step
    mu_nom <- matrix(0, nrow = cluster_num, ncol = dim_data)
    mu_denom<-array(0, dim = c(dim_data, dim_data, cluster_num))
    
    for (g in seq_along(g_clusternum)) {
        for (k in seq_len(cluster_num)) {
            cov_inv <- solve(covh[, , k])
            mu_nom[k, ] <- mu_nom[k, ] + zi[g, k] * mg[g] * (mug[g, ] %*% cov_inv)
            mu_denom[,,k] <- mu_denom[,,k] + zi[g, k] * mg[g] * cov_inv
        }
    }
    
    primary_seq<-seq(1,(log2(cluster_num)+1))
    mus<-muh[seq_len(log2(cluster_num)+1),]
    mus_nonprimary<-muh[(log2(cluster_num)+2):nrow(muh),]
    
    nonprimary_seq<-seq((log2(cluster_num)+2),cluster_num)
    
    for (i in primary_seq){
        weights_other <- weights[,-i]%*%mus[-i,]-mus_nonprimary
        mu_nom[i,] <- mu_nom[i,] - (2*lambdas[-1]*t(weights[,i]))%*% weights_other
        mu_denom[,,i] <- mu_denom[,,i] + ((2*lambdas[-1]*t(weights[,i]))%*%weights[,i])[1,1]*diag(1,nrow=dim_data)
    }
    
    mu_nom[1,] <- mu_nom[1,] + 2 * lambdas[1] * neg_assum
    mu_denom[,,1] <- mu_denom[,,1] + 2 * lambdas[1] * diag(1, nrow = dim_data) 
    
    mu_nom[-seq_len(log2(cluster_num)+1),]<- mu_nom[-seq_len(log2(cluster_num)+1),] + 2*lambdas[-1]*(weights%*%mus)
    
    diag_tmp<-array(0, dim = c(dim_data, dim_data, length(nonprimary_seq)))
    for (j in seq_along(nonprimary_seq)) {
        diag_tmp[,,j] <- diag(dim_data)  
    }
    
    if (cluster_num==4){
        mu_denom[,,-seq_len(log2(cluster_num)+1)]<- mu_denom[,,-seq_len(log2(cluster_num)+1)] + (2*lambdas[-1]*diag(2))
    }else{
        mu_denom[,,-seq_len(log2(cluster_num)+1)]<- mu_denom[,,-seq_len(log2(cluster_num)+1)] + (2*lambdas[-1]*diag_tmp)
    }
    
    for (h in seq_len(cluster_num)){
        muh[h,] <- mu_nom[h,]%*%solve(mu_denom[,,h])
    }
    
    return(muh)
}



#' Internal Function 9
#'
#' This function calculates mu in M-step of EM algorithm
#' @param zi the expected log-likelihood found on the E step
#' @param g_clusternum cluster labels from base clustering 
#' @param cluster_num The expected maximum number of clusters
#' @param dim_data the dimension of the dataset
#' @param muh the matrix of cluster centers at level l
#' @param mg cluster sizes of base clustering result
#' @param mug the matrix of cluster centers at level l+1
#' @param covg the covariance matrix of clusters at level l+1
#' @return covh the covariance matrix of clusters at level l in the EM algorithm
#' @keywords internal
mstep_cov<-function(cluster_num,dim_data,g_clusternum,zi,mg,covg,mug,muh){
    # Update covariances
    sigma_nom <- matrix(0, nrow = dim_data, ncol = dim_data)
    sigma_denom <- 0
    covh <- array(cov(mug), dim = c(dim_data, dim_data, cluster_num))
    for (k in seq_len(cluster_num)) {
        for (g in seq_len(length(g_clusternum))) {
            sigma_nom <- sigma_nom + zi[g, k] * mg[g] * covg[, , g] + zi[g, k] * mg[g] * (mug[g, ] - muh[k, ]) %*% t(mug[g, ] - muh[k, ])
            sigma_denom <- sigma_denom + zi[g, k] * mg[g]
        }
        covh[, , k] <- sigma_nom / sigma_denom
    }
    return(covh)
}

#' Internal Function 10
#'
#' This function merges the excess clusters given by the base clustering
#' @param data A matrix or data frame of fluorescence intensities in each channel. Each row represents each partitions, and each column each channel.
#' @param cluster_num The expected maximum number of clusters
#' @param base_cluster base clustering results before merging
#' @param eps the convergence threshold
#' @param max_iter maximum number of iterations
#' @param lambdas The penalty terms for the deviation from the expected cluster centers. Higher \code{lambdas} penalizes the deviation more.
#' @param coefs The coefficients to adjust for the expected cluster centers. The default is 1 which can be used for common assay designs and has
#' to be modified for special assays such as competing assays.
#' @return A list of membership probability, cluster center, merging probability
#' @keywords internal
HMM_merge <- function(data, cluster_num, base_clust, eps = 10^(-10), max_iter = 1000, lambdas = rep(2, 2), coefs = rep(1, 2)) {
    init_results<-GMM_init(data, cluster_num, base_clust, coefs)
    pih<-init_results[[1]]
    muh<-init_results[[2]]
    covh<-init_results[[3]]
    g_clusternum<-init_results[[4]]
    mg<-init_results[[5]]
    mug<-init_results[[6]]
    covg<-init_results[[7]]
    weights<-init_results[[8]]
    neg_assum<-init_results[[9]]
    
    dim_data <- ncol(data)
    change_ests <- NULL
    
    # Start the EM algorithm
    for (j in seq_len(max_iter)) {
        print(j)
        if (j >= max_iter) {
            message("Warning: the algorithm fails to converge")
        }
        
        zi<-estep(g_clusternum,cluster_num,pih,muh,covh,mg,mug,covg)
        
        pih <- apply(zi,2,sum)/length(g_clusternum)
        pih[which(pih==0)]<-(10^(-10))
        
        muh<-mstep_mu(zi,g_clusternum,dim_data,cluster_num,weights,muh,covh,mg,mug,neg_assum,lambdas,coefs)
        
        covh<-mstep_cov(cluster_num,dim_data,g_clusternum,zi,mg,covg,mug,muh)
        
        
        # Change the likelihood
        change_ests <- c(change_ests, sum(apply(zi, 1, function(x) log(sum(x)))))
        
        # Stopping criteria
        if (j > 1) {
            if (sum(abs(muh - old_mu)) < eps) break
        }
        
        old_mu <- muh
    }
    
    return(list(zi, muh, pih))
}


#' Internal Function 11
#'
#' This function outputs all combinations of primary targets
#' @param cluster_num The expected maximum number of clusters
#' @return A matrix of all combinations of primary targets
#' @keywords internal

cluster_selection<-function(cluster_num){
    primary_tar<-log2(cluster_num)
    all_combinations <- list()
    mat_coef<-NULL
    
    comb_seq<-seq(1,primary_tar)
    coef_tmp<-rep(0,primary_tar)
    
    for (i in comb_seq){
        comb_tmp<-combn(seq_len(primary_tar), i)
        all_combinations[[i]] <- comb_tmp
        comb_ncol<-ncol(comb_tmp)
        for (j in seq_len(comb_ncol)){
            mat_coef_tmp<-coef_tmp
            mat_coef_tmp[comb_tmp[,j]]<-1
            mat_coef<-rbind(mat_coef,mat_coef_tmp)
        }
    }
    
    mat_coef<-rbind(rep(0,primary_tar),mat_coef)
    return(mat_coef)
}
