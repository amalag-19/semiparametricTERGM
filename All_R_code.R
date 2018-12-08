########################################################################################################
## Contains the most updated working code for semiparametric VEM and SEM, undirected and directed networks for all K.
#########################################################################################################
## Loading the required packages
require(lda)
library(quadprog)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(ggplot2)
library(abind)
library(reshape)
library(combinat)
require(network)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#########################################################################################################
## Defining a wrapper function for EM undirected case to run for different number of clusters and bandwidths
wrapper_EM_undir<-function(sim.net,nclust,bandwidth,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(EMundir)
  library(EMundirK1)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-matrix(NA_real_,T_grid,K)
    for (t in 1:T_grid){
      kernel_vec<-rep(0,T_data)
      for (temp_index in 1:T_data){
        kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      nonzero_ids<-which(kernel_vec!=0)
      nonzero_ids_len<-length(nonzero_ids)
      grad_vec_kern<-matrix(0,T_data,K)
      if (nonzero_ids_len!=0){
        grad_vec<-grad_EM_undir(theta_u=as.vector(theta.curr[t,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len)
        for (temp_index in 1:nonzero_ids_len){
          grad_vec_kern[nonzero_ids[temp_index],]<-grad_vec[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
        }
      }
      gradient<-colSums(grad_vec_kern)
      hess<-hess_EM_undir(theta_u=as.vector(theta.curr[t,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      if(sum(hess)==0){
        theta.next[t,]<-as.vector(theta.curr[t,])
      } else{theta.next[t,]<-as.vector(theta.curr[t,])-as.vector(solve(hess)%*%gradient)}
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the indices for actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(T_grid,K,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta matrix
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_undir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Defining Model Selection functions based on converged paramters for K>=2.
  
  ## Defining a function to calculate the estimate of H_K
  H_K<-function(gamma,theta,N,K,T_grid){
    H_K_mat<-matrix(NA_real_,K,K)
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    ## Filling the upper triangular matrix of H_K since it is symmetric
    for (k in 1:K){
      for (l in k:K){
        if (k==l){
          t1<-0
          for (t  in 1:T_grid){
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                indicator_1<-gamma[i,k]>(1/K)
                indicator_2<-gamma[j,k]>(1/K)
                if(indicator_1|indicator_2){
                  indicator<-indicator_1+indicator_2
                  exp_val<-exp(2*theta_data[t,k])
                  t1<-t1+((exp_val/((1+exp_val)^2))*(indicator^2))
                }
              }
            }
          }
          H_K_mat[k,k]<-t1
        }else if (k!=l){
          t1<-0
          for (t  in 1:T_grid){
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                if(((gamma[i,k]>(1/K))&(gamma[j,l]>(1/K)))|((gamma[i,l]>(1/K))&(gamma[j,k]>(1/K)))){
                  exp_val<-exp(theta_data[t,k]+theta_data[t,l])
                  t1<-t1+(exp_val/((1+exp_val)^2))
                }
              }
            }
          }
          H_K_mat[k,l]<-t1
          H_K_mat[l,k]<-t1
        }
      }
    }
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K<-function(gamma,theta,network,N,K,T_grid,grid_ids){
    V_K_mat<-matrix(0,K,K)
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    for (t in 1:T_grid){
      u_vec<-rep(0,K)
      for (k in 1:K){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            indicator_i<-(gamma[i,k]>(1/K))
            indicator_j<-(gamma[j,k]>(1/K))
            if(indicator_i|indicator_j){
              if(indicator_i&indicator_j){
                exp_val<-exp(2*theta_data[t,k])
                u_vec[k]<-u_vec[k]+((network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))*(indicator_i+indicator_j))
              }else if(indicator_i&(!indicator_j)){
                clusterid_j<-which.max(gamma[j,])
                exp_val<-exp(theta_data[t,k]+theta_data[t,clusterid_j])
                u_vec[k]<-u_vec[k]+(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))
              }else if(indicator_j&(!indicator_i)){
                clusterid_i<-which.max(gamma[i,])
                exp_val<-exp(theta_data[t,k]+theta_data[t,clusterid_i])
                u_vec[k]<-u_vec[k]+(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))
              }
            }
          }
        }
      }
      V_K_mat<-V_K_mat+u_vec%*%t(u_vec)
    }
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(gamma,theta,network,N,K,T_grid,grid_ids){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-which.max(gamma[i,])
          cluster_id_j<-which.max(gamma[j,])
          exp_val<-exp(theta_data[t,cluster_id_i]+theta_data[t,cluster_id_j])
          cl_val<-cl_val+((network[i,j,grid_ids[t]]*(theta_data[t,cluster_id_i]+theta_data[t,cluster_id_j]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC<-function(gamma,theta,network,N,K,T_data,T_grid,bandwidth,grid_ids){
    H_K_mat<-H_K(gamma = gamma,theta = theta,N = N,K = K,T_grid = T_grid)
    V_K_mat<-V_K(gamma = gamma,theta = theta,network = network, N = N,K = K,T_grid = T_grid,grid_ids=grid_ids)
    d_K<-sum(diag(solve(H_K_mat+(10^(-10)*diag(K)))%*%(V_K_mat)))
    #d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    t1<-(-2*(cond_loglik(gamma = gamma,theta = theta,network = network,N = N,K = K,T_grid = T_grid,grid_ids=grid_ids)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-rep(NA_real_,T_grid)
    for (t in 1:T_grid){
      grad_vec<-grad_EM_undir_K1(theta_u=theta.curr[t], network=network, N=N, K=K, T_data=T_data)
      grad_vec_kern<-rep(NA_real_,T_data)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-sum(grad_vec_kern)
      hess<-hess_EM_undir_K1(theta_u=theta.curr[t], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    theta<-matrix(NA_real_,T_grid,n_iter)
    theta[,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_undir_K1(theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  ## Defining Model Selection functions for K=1 based on converged paramters.
  
  ## Defining a function to calculate the estimate of H_K
  H_K1<-function(theta,N,T_grid){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    H_K_val<-0
    for (t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          H_K_val<-H_K_val+(4*(exp_val/((1+exp_val)^2)))
        }
      }
    }
    return(H_K_val)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K1<-function(theta,network,N,T_grid,grid_ids){
    V_K_val<-0
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    for (t in 1:T_grid){
      u_val<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          u_val<-u_val+(2*(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val))))
        }
      }
      V_K_val<-V_K_val+(u_val^2)
    }
    return(V_K_val)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC_K1<-function(theta,network,N,K,T_data,T_grid,bandwidth,grid_ids){
    H_K_val<-H_K1(theta = theta,N = N,T_grid = T_grid)
    V_K_val<-V_K1(theta = theta,network = network, N = N,T_grid = T_grid,grid_ids=grid_ids)
    d_K<-V_K_val/H_K_val
    t1<-(-2*(cond_loglik_K1(theta = theta,network = network,N = N,T_grid = T_grid,grid_ids=grid_ids)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE function
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  ## Defining the actual data points and grid points
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  K<-nclust ## Defining the number of clusters
  
  #################################################
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
    for(i in 1:N){
      gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
      gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-rep(0,T_grid)
    #debug(iterator_K1)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K)#pi.update(gamma.curr = gamma.start,N=N,K=K)
    start[[3]]<-matrix(0,T_grid,K)
    #debug(iterator)
    #ptm<-proc.time()
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
    #print(proc.time()-ptm)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,n_last]
    BIC_val<-BIC_K1(theta = param_converge,network = sim.net,N = N,T_data = T_data,T_grid = T_grid,bandwidth = h,grid_ids=grid_ids)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = rep(1,N),cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_theta)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    BIC_val<-BIC(gamma = param_converge[[1]], theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data,T_grid = T_grid,bandwidth = h,grid_ids=grid_ids)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][,K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id<-which.min(RASE_theta_vec)
        RASE_theta<-RASE_theta_vec[permute_true_id]
        output_list<-list(param_converge,BIC_val,RI_val,RASE_theta)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper cross validation function to run for fixed number of clusters and a given bandwidth
wrapper_EM_undir_CV<-function(sim.net,nclust,bandwidth,fold_list,folds){
  ## This is the combined final updated and working code for NCV for EM directed case for all K
  #########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(EMundir)
  library(EMundirK1)
  library(Matrix)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids,test_node_set_len){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    quad_lin_coeff<-gamma_update_undir_CV(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
    
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of the theta with dimensions T_grid*K*2
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
    theta.next<-array(NA_real_,dim=c(T_grid,K))
    for (t in 1:T_grid){
      kernel_vec<-rep(0,T_data)
      for (temp_index in 1:T_data){
        kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      nonzero_ids<-which(kernel_vec!=0)
      nonzero_ids_len<-length(nonzero_ids)
      grad_vec_kern<-matrix(0,T_data,K)
      if (nonzero_ids_len!=0){
        grad_vec<-grad_EM_undir_CV(theta_u=as.matrix(theta.curr[t,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,test_node_set_len = test_node_set_len)
        for (temp_index in 1:nonzero_ids_len){
          grad_vec_kern[nonzero_ids[temp_index],]<-grad_vec[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
        }
      }
      gradient<-as.vector(colSums(grad_vec_kern))
      hess<-hess_EM_undir_CV(theta_u=as.matrix(theta.curr[t,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
     
      #theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)
      if(sum(hess)==0){
        theta.next[t,]<-theta.curr[t,]
      } else{theta.next[t,]<-theta.curr[t,]-as.vector(solve(hess)%*%gradient)}
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the indices for actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(T_grid,K,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)&(iter_index<25)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta array
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len = test_node_set_len)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_undir_CV(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  ##########################################################################################################################################################################################################
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(theta,network,N,K,T_data,cluster_ids,test_node_set_len){
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    cl_val<-0
    for(t in 1:T_data){
      for(i in 1:(test_node_set_len-1)){
        for(j in (i+1):test_node_set_len){
          exp_val<-exp(theta[t,cluster_ids[i]]+theta[t,cluster_ids[j]])
          cl_val<-cl_val+((network[i,j,t]*(theta[t,cluster_ids[i]]+theta[t,cluster_ids[j]]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining the functions for K=1
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
    theta.next<-rep(NA_real_,T_grid)
    for (t in 1:T_grid){
      grad_vec<-grad_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data,test_node_set_len=test_node_set_len)
      grad_vec_kern<-rep(NA_real_,T_data)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-sum(grad_vec_kern)
      hess<-hess_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
      theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-matrix(NA_real_,N,n_iter)
    pi<-rep(NA_real_,n_iter)
    theta<-matrix(NA_real_,T_grid,n_iter)
    gamma[,1]<-start[[1]]
    pi[1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      #ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,iter_index]<-rep(1,N)
      ## Updating the pi vector
      pi[iter_index]<-1
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], pi=pi[iter_index],gamma=gamma[,iter_index], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len=test_node_set_len) 
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_K1(gamma=gamma[,iter_index], pi=pi[iter_index], theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      #print(iter_index)
      #print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  ##########################################################################################################################################################################################################
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids,test_node_set_len){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:test_node_set_len){
        if(i!=test_node_set_len){
          for(j in (i+1):test_node_set_len){
            exp_val<-exp(2*theta_data[t])
            cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
          }
        } else if(i==test_node_set_len){
          exp_val<-exp(2*theta_data[t])
          cl_val<-cl_val+((network[i,test_node_set_len,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  #################################################
  ## Defining the parameters 
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  #T_grid<-T_data
  
  ## Defining the actual data points and grid points 
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  K<-nclust ## Defining the number of clusters
  
  ## Initializing the cond_loglik_val vector to store the loss corresponding to each fold
  cond_loglik_val<-rep(NA_real_,folds)
  ## Loop over different folds for cross validation to choose bandwidth
  for (fold_index in 1:folds){
    test_node_set<-sort(fold_list[[fold_index]])
    train_list<-fold_list
    train_list[[fold_index]]<-NULL
    train_node_set<-sort(unlist(train_list))
    permute_set<-c(test_node_set,train_node_set)
    sim.net_permute<-array(NA_real_,dim = dim(sim.net))
    for (i in 1:T_data){
      sim.net_permute[,,i]<-as(as.integer(permute_set), "pMatrix")%*%sim.net[,,i]%*%t(as(as.integer(permute_set), "pMatrix"))
    }
    
    #################################################
    ## Setting the initial value of gamma
    if(K==1){
      gamma.start <- rep(1,N)
    }else{
      ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
      set.seed((2))
      MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
      gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
      ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
      for(i in 1:N){
        gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
        gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
      }
    }
    
    #################################################
    ## Defining the starting values of the iterator and running the main algorithm
    if(K==1){
      start<-matrix(0,T_grid,2)
      param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
    }else{
      start<-list()
      start[[1]]<-gamma.start
      start[[2]]<-rep(1/K,K)
      start[[3]]<-array(0,dim=c(T_grid,K))
      #debug(iterator)
      param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids,test_node_set_len = length(test_node_set))
    }
    
    #################################################
    ## extracting the coverged parameter values and calculating BIC
    n_iter=1
    indicator_last<-0
    while(indicator_last==0){
      if(K==1){
        temp<-is.na(param[1,1,n_iter])
      }else{temp<-is.na(param[[1]][1,1,n_iter])}
      if(temp==TRUE){
        n_last<-n_iter-1
        indicator_last<-1
      }
      n_iter<-n_iter+1
    }
    param_converge<-list()
    if(K==1){
      param_converge[[1]]<-param[[1]][,n_last]
      param_converge[[2]]<-param[[2]][n_last]
      param_converge[[3]]<-param[[3]][,n_last]
      cond_loglik_val[fold_index]<-cond_loglik_K1(theta = param_converge[[3]],network = sim.net_permute,N = N,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = length(test_node_set))
    }else{
      param_converge[[1]]<-param[[1]][,,n_last]
      param_converge[[2]]<-param[[2]][,n_last]
      param_converge[[3]]<-param[[3]][,,n_last]
      cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        cluster_id<-which.max(param_converge[[1]][x,])
        return(cluster_id)
      }))
      cond_loglik_val[fold_index]<-cond_loglik(theta = param_converge[[3]],network = sim.net_permute,N = N,K = K,T_data = T_data,cluster_ids=cluster_ids,test_node_set_len = length(test_node_set))
    }
  }
  cond_loglik_sum<--sum(cond_loglik_val)
  return(cond_loglik_sum)
}

#########################################################################################################
## Defining a function band_wrapper for doing cross validation and choosing bandwidth for a given repeat
band_wrapper_EM_undir<-function(sim.net,bandwidths,K,folds){
  ###################################################
  ## Defining a wrapper cross validation function to run for fixed number of clusters and a given bandwidth
  wrapper_EM_undir_CV<-function(sim.net,nclust,bandwidth,fold_list,folds){
    ## This is the combined final updated and working code for NCV for EM directed case for all K
    #########################################################################################################
    ## Loading the required packages
    require(lda)
    library(quadprog)
    library(combinat)
    library(EMundir)
    library(EMundirK1)
    library(Matrix)
    
    ######################################################################################################
    ## Defining the epanechnikov kernel function
    epan<-function(y){
      return(0.75*(1-y^2)*(abs(y)<=1))
    }
    
    #################################################
    ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
    ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
    ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
    gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids,test_node_set_len){
      gamma.next<-matrix(NA_real_,N,K)
      constraint_matrix<-matrix(NA_real_,2+K,K)
      constraint_matrix[1,]<-rep(1,K)
      constraint_matrix[2,]<-rep(-1,K)
      constraint_matrix[3:(K+2),]<-diag(K)
      constraint_vector<-c(1,-1,rep(0,K))
      quad_lin_coeff<-gamma_update_undir_CV(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
      
      for (i in 1:N){
        gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
      }
      
      ## normalizing gamma_i. deal with case outside (0,1) later
      gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
        x_norm<-x/sum(x)
        return(x_norm)
      }))
      return(gamma.next)
    }
    
    #################################################
    ## Defining a function to update K*1 vector of pi
    pi.update<-function(gamma.curr,N,K){
      pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
      ## normalization of pi
      pi.next<-pi.next/sum(pi.next)
      return(pi.next)
    }
    
    #################################################
    ## Defining the update function of the theta with dimensions T_grid*K*2
    theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
      theta.next<-array(NA_real_,dim=c(T_grid,K))
      for (t in 1:T_grid){
        kernel_vec<-rep(0,T_data)
        for (temp_index in 1:T_data){
          kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
        }
        nonzero_ids<-which(kernel_vec!=0)
        nonzero_ids_len<-length(nonzero_ids)
        grad_vec_kern<-matrix(0,T_data,K)
        if (nonzero_ids_len!=0){
          grad_vec<-grad_EM_undir_CV(theta_u=as.matrix(theta.curr[t,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,test_node_set_len = test_node_set_len)
          for (temp_index in 1:nonzero_ids_len){
            grad_vec_kern[nonzero_ids[temp_index],]<-grad_vec[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
          }
        }
        gradient<-as.vector(colSums(grad_vec_kern))
        hess<-hess_EM_undir_CV(theta_u=as.matrix(theta.curr[t,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
        
        #theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)
        if(sum(hess)==0){
          theta.next[t,]<-theta.curr[t,]
        } else{theta.next[t,]<-theta.curr[t,]-as.vector(solve(hess)%*%gradient)}
      }
      return(theta.next)
    }
    
    #################################################
    ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
    iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
      N<-dim(network)[1]
      T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
      ## adjacency matrix
      T_grid<-length(grid_ids)  ## Number of grid points
      
      ## Defining the indices for actual data points
      data_ids<-1:T_data
      
      ## initializing the arrays for parameters
      gamma<-array(NA_real_,dim=c(N,K,n_iter))
      pi<-matrix(NA_real_,K,n_iter)
      theta<-array(NA_real_,dim=c(T_grid,K,n_iter))
      gamma[,,1]<-start[[1]]
      pi[,1]<-start[[2]]
      theta[,,1]<-start[[3]]
      ## iterations
      ## Defining the iteration index
      iter_index<-2
      ## Initializing the error
      error<-Inf
      ## Initializing the current ELBO values over the whole grid
      ELBO_grid.curr<-rep(10^10,T_grid)
      while((error>thres)&(iter_index<25)){
        ## Starting the stopwatch to calculate the per iteration time
        ptm<-proc.time()
        ## Updating the N*K gamma matrix i.e. variational variational parameters
        gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
        ## Updating the pi vector
        pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
        ## Updating the theta array
        theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len = test_node_set_len)
        ELBO_grid.prev<-ELBO_grid.curr
        ELBO_grid.curr<-rep(NA_real_,T_grid)
        for(t in 1:T_grid){
          ELBO_grid.curr[t]<-ELBO_conv_EM_undir_CV(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
        }
        error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
        print(error)
        print(iter_index)
        print(proc.time()-ptm)
        iter_index<-iter_index+1
      }
      return(list(gamma,pi,theta))
    }
    
    ##########################################################################################################################################################################################################
    ## Defining a function to calculate the estimate of conditional log likelihood
    cond_loglik<-function(theta,network,N,K,T_data,cluster_ids,test_node_set_len){
      ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
      cl_val<-0
      for(t in 1:T_data){
        for(i in 1:(test_node_set_len-1)){
          for(j in (i+1):test_node_set_len){
            exp_val<-exp(theta[t,cluster_ids[i]]+theta[t,cluster_ids[j]])
            cl_val<-cl_val+((network[i,j,t]*(theta[t,cluster_ids[i]]+theta[t,cluster_ids[j]]))-(log(1+exp_val)))
          }
        }
      }
      return(cl_val)
    }
    
    #########################################################################################################
    #########################################################################################################
    ## Defining the functions for K=1
    #################################################
    ## Defining the update function of full matrix theta with dimensions T_grid*K
    theta.update_K1<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
      theta.next<-rep(NA_real_,T_grid)
      for (t in 1:T_grid){
        grad_vec<-grad_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data,test_node_set_len=test_node_set_len)
        grad_vec_kern<-rep(NA_real_,T_data)
        for (temp_index in 1:T_data){
          grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
        }
        gradient<-sum(grad_vec_kern)
        hess<-hess_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
        theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
      }
      return(theta.next)
    }
    
    #################################################
    ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
    iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
      N<-dim(network)[1]
      T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
      ## adjacency matrix
      T_grid<-length(grid_ids)  ## Number of grid points
      
      ## Defining the actual data points
      data_ids<-1:T_data
      
      ## initializing the arrays for parameters
      gamma<-matrix(NA_real_,N,n_iter)
      pi<-rep(NA_real_,n_iter)
      theta<-matrix(NA_real_,T_grid,n_iter)
      gamma[,1]<-start[[1]]
      pi[1]<-start[[2]]
      theta[,1]<-start[[3]]
      ## iterations
      ## Defining the iteration index
      iter_index<-2
      ## Initializing the error
      error<-Inf
      ## Initializing the current ELBO values over the whole grid
      ELBO_grid.curr<-rep(10^10,T_grid)
      while(error>thres){
        ## Starting the stopwatch to calculate the per iteration time
        #ptm<-proc.time()
        ## Updating the N*K gamma matrix i.e. variational variational parameters
        gamma[,iter_index]<-rep(1,N)
        ## Updating the pi vector
        pi[iter_index]<-1
        ## Updating the theta matrix
        theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], pi=pi[iter_index],gamma=gamma[,iter_index], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len=test_node_set_len) 
        ELBO_grid.prev<-ELBO_grid.curr
        ELBO_grid.curr<-rep(NA_real_,T_grid)
        for(t in 1:T_grid){
          ELBO_grid.curr[t]<-ELBO_conv_K1(gamma=gamma[,iter_index], pi=pi[iter_index], theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
        }
        error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
        #print(iter_index)
        #print(proc.time()-ptm)
        iter_index<-iter_index+1
      }
      return(list(gamma,pi,theta))
    }
    
    ##########################################################################################################################################################################################################
    ## Defining a function to calculate the estimate of conditional log likelihood
    cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids,test_node_set_len){
      ## Assuming grid spacing is half the data spacing.
      theta_data<-theta
      cl_val<-0
      for(t in 1:T_grid){
        for(i in 1:test_node_set_len){
          if(i!=test_node_set_len){
            for(j in (i+1):test_node_set_len){
              exp_val<-exp(2*theta_data[t])
              cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
            }
          } else if(i==test_node_set_len){
            exp_val<-exp(2*theta_data[t])
            cl_val<-cl_val+((network[i,test_node_set_len,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
          }
        }
      }
      return(cl_val)
    }
    
    #################################################
    ## Defining the parameters 
    N<-dim(sim.net)[1] ## Number of nodes from the network
    ## Total time points in the network for which we have data in the form of adjacency matrix
    T_data<-dim(sim.net)[3]
    
    #T_grid<-T_data
    
    ## Defining the actual data points and grid points 
    data_ids<-(1:T_data)
    #grid_ids<-1:T_data
    grid_ids<-seq(1,T_data,by = 1)
    
    ## Defining the total number of grid points
    T_grid<-length(grid_ids)
    
    # Defining the epanechnikov kernel bandwidth: h
    h<-bandwidth
    
    K<-nclust ## Defining the number of clusters
    
    ## Initializing the cond_loglik_val vector to store the loss corresponding to each fold
    cond_loglik_val<-rep(NA_real_,folds)
    ## Loop over different folds for cross validation to choose bandwidth
    for (fold_index in 1:folds){
      test_node_set<-sort(fold_list[[fold_index]])
      train_list<-fold_list
      train_list[[fold_index]]<-NULL
      train_node_set<-sort(unlist(train_list))
      permute_set<-c(test_node_set,train_node_set)
      sim.net_permute<-array(NA_real_,dim = dim(sim.net))
      for (i in 1:T_data){
        sim.net_permute[,,i]<-as(as.integer(permute_set), "pMatrix")%*%sim.net[,,i]%*%t(as(as.integer(permute_set), "pMatrix"))
      }
      
      #################################################
      ## Setting the initial value of gamma
      if(K==1){
        gamma.start <- rep(1,N)
      }else{
        ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
        set.seed((2))
        MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
        gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
        ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
        for(i in 1:N){
          gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
          gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
        }
      }
      
      #################################################
      ## Defining the starting values of the iterator and running the main algorithm
      if(K==1){
        start<-matrix(0,T_grid,2)
        param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
      }else{
        start<-list()
        start[[1]]<-gamma.start
        start[[2]]<-rep(1/K,K)
        start[[3]]<-array(0,dim=c(T_grid,K))
        #debug(iterator)
        param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids,test_node_set_len = length(test_node_set))
      }
      
      #################################################
      ## extracting the coverged parameter values and calculating BIC
      n_iter=1
      indicator_last<-0
      while(indicator_last==0){
        if(K==1){
          temp<-is.na(param[1,1,n_iter])
        }else{temp<-is.na(param[[1]][1,1,n_iter])}
        if(temp==TRUE){
          n_last<-n_iter-1
          indicator_last<-1
        }
        n_iter<-n_iter+1
      }
      param_converge<-list()
      if(K==1){
        param_converge[[1]]<-param[[1]][,n_last]
        param_converge[[2]]<-param[[2]][n_last]
        param_converge[[3]]<-param[[3]][,n_last]
        cond_loglik_val[fold_index]<-cond_loglik_K1(theta = param_converge[[3]],network = sim.net_permute,N = N,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = length(test_node_set))
      }else{
        param_converge[[1]]<-param[[1]][,,n_last]
        param_converge[[2]]<-param[[2]][,n_last]
        param_converge[[3]]<-param[[3]][,,n_last]
        cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
          cluster_id<-which.max(param_converge[[1]][x,])
          return(cluster_id)
        }))
        cond_loglik_val[fold_index]<-cond_loglik(theta = param_converge[[3]],network = sim.net_permute,N = N,K = K,T_data = T_data,cluster_ids=cluster_ids,test_node_set_len = length(test_node_set))
      }
    }
    cond_loglik_sum<--sum(cond_loglik_val)
    return(cond_loglik_sum)
  }
  
  ###################################################
  ## Defining a function to get N_Vfolds used as a input to wrapper_CV type of functions while doing NCV
  get_fold_list<-function(sim.net,folds){
    ## Dividing the nodes set into V folds equal subsets where V = 3 is most commonly used
    ## randomly sampling the nodes set
    N<-dim(sim.net)[1] ## Number of nodes from the network
    N_set_rand<-sample(N)
    N_Vfolds<-list()
    step_size<-floor(N/folds)
    for (i in 1:folds){
      N_Vfolds[[i]]<-N_set_rand[(1+(i-1)*step_size):(i*step_size)]
    }
    if((N%%folds)!=0){
      for (i in 1:(N%%folds)){
        N_Vfolds[[i]]<-c(N_Vfolds[[i]],N_set_rand[i+(step_size*folds)])
      }
    }
    return(N_Vfolds)
  }
  
  ###################################################
  ## Initializing the conditional log-likelihood vector for different bandwidths
  cond_loglik_vec<-rep(NA_real_,length(bandwidths))
  for (band_index in 1:length(bandwidths)){
    cond_loglik_vec[band_index]<- wrapper_EM_undir_CV(sim.net = sim.net,nclust = K,bandwidth = bandwidths[band_index],fold_list = get_fold_list(sim.net = sim.net,folds = folds),folds = folds)
  }
  bandwidth_chosen<-bandwidths[which.min(cond_loglik_vec)]
  return(bandwidth_chosen)
}

#########################################################################################################
## Defining a wrapper function for EM directed case to run for different number of clusters and bandwidths
wrapper_EM_dir<-function(sim.net,nclust,bandwidth,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM directed case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(EMdir)
  library(EMdirK1)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    quad_lin_coeff<-gamma_update_dir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids)
    
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of the theta with dimensions T_grid*K*2
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-array(NA_real_,dim=c(T_grid,K,2))
    for (t in 1:T_grid){
      kernel_vec<-rep(0,T_data)
      for (temp_index in 1:T_data){
        kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      nonzero_ids<-which(kernel_vec!=0)
      nonzero_ids_len<-length(nonzero_ids)
      grad_vec_kern_oe<-matrix(0,T_data,K)
      grad_vec_kern_re<-matrix(0,T_data,K)
      if (nonzero_ids_len!=0){
        grad_oe<-grad_EM_dir_oe(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len)
        grad_re<-grad_EM_dir_re(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len)
        for (temp_index in 1:nonzero_ids_len){
          grad_vec_kern_oe[nonzero_ids[temp_index],]<-grad_oe[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
          grad_vec_kern_re[nonzero_ids[temp_index],]<-grad_re[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
        }
      }
      gradient_oe<-as.vector(colSums(grad_vec_kern_oe))
      gradient_re<-as.vector(colSums(grad_vec_kern_re))
      gradient<-c(gradient_oe,gradient_re)
      hess_oe<-hess_EM_dir_oe(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess_re<-hess_EM_dir_re(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess_oe_re<-hess_EM_dir_oe_re(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess<-matrix(NA_real_,2*K,2*K)
      hess[1:K,1:K]<-hess_oe
      hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
      hess[(1:K),((K+1):(2*K))]<-hess_oe_re
      hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
      #theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)
      if(sum(hess)==0){
        theta.next[t,,]<-theta.curr[t,,]
      } else{theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)}
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the indices for actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(T_grid,K,2,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)&(iter_index<25)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta array
      theta[,,,iter_index]<-theta.update(theta.curr=theta[,,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_dir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Model Selection functions for K>=2
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(theta,network,N,K,T_data,cluster_ids){
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    cl_val<-0
    for(t in 1:T_data){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val_1<-exp(theta[t,cluster_ids[i],1])
          exp_val_2<-exp(theta[t,cluster_ids[j],1])
          exp_val_3<-exp(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2])
          indicator_10<-(network[i,j,t]==1)&(network[j,i,t]==0)
          indicator_01<-(network[i,j,t]==0)&(network[j,i,t]==1)
          indicator_11<-(network[i,j,t]==1)&(network[j,i,t]==1)
          cl_val<-cl_val+((indicator_10*theta[t,cluster_ids[i],1])+(indicator_01*theta[t,cluster_ids[j],1])+(indicator_11*(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2]))-(log(1+exp_val_1+exp_val_2+exp_val_3)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of H_K
  H_K<-function(theta,N,K,T_data,cluster_ids){
    H_K_mat<-matrix(NA_real_,2*K,2*K)
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    H_K_oe<-HK_dir_oe(theta=theta, N=N, K=K, T_data=T_data, cluster_ids=cluster_ids)
    H_K_re<-HK_dir_re(theta=theta, N=N, K=K, T_data=T_data, cluster_ids=cluster_ids)
    H_K_oe_re<-HK_dir_oe_re(theta=theta, N=N, K=K, T_data=T_data, cluster_ids=cluster_ids)
    H_K_mat[1:K,1:K]<-H_K_oe
    H_K_mat[((K+1):(2*K)),((K+1):(2*K))]<-H_K_re
    H_K_mat[(1:K),((K+1):(2*K))]<-H_K_oe_re
    H_K_mat[((K+1):(2*K)),(1:K)]<-t(H_K_oe_re)
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K<-function(theta,network,N,K,T_data,cluster_ids){
    V_K_mat<-matrix(0,2*K,2*K)
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    u_mat_oe<-VK_dir_oe(theta=theta, network=network, N=N, K=K, T_data=T_data,cluster_ids=cluster_ids)
    u_mat_re<-VK_dir_re(theta=theta, network=network, N=N, K=K, T_data=T_data,cluster_ids=cluster_ids)
    for (t in 1:T_data){
      u_vec<-c(as.vector(u_mat_oe[t,]),as.vector(u_mat_re[t,]))
      V_K_mat<-V_K_mat+u_vec%*%t(u_vec)
    }
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC<-function(theta,network,N,K,T_data,bandwidth,cluster_ids){
    H_K_mat<-H_K(theta = theta,N = N,K = K,T_data = T_data,cluster_ids = cluster_ids)
    V_K_mat<-V_K(theta = theta,network = network, N = N,K = K,T_data = T_data,cluster_ids=cluster_ids)
    d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    cond_loglik_val<-cond_loglik(theta = theta,network = network,N = N,K = K,T_data = T_data,cluster_ids=cluster_ids)
    t1<-(-2*cond_loglik_val)
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(cond_loglik_val,t2,BIC_val))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of theta_oe/ie with dimensions T_grid
  theta.update_K1<-function(theta.curr,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-matrix(NA_real_,T_grid,2)
    for (t in 1:T_grid){
      grad_oe_ie<-grad_EM_dir_K1_oe(theta_u=theta.curr[t,], network=network, N=N, K=K, T_data=T_data)
      grad_re<-grad_EM_dir_K1_re(theta_u=theta.curr[t,], network=network, N=N, K=K, T_data=T_data)
      grad_vec_kern<-matrix(NA_real_,T_data,2)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index,1]<-grad_oe_ie[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
        grad_vec_kern[temp_index,2]<-grad_re[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-as.vector(colSums(grad_vec_kern))
      hess_oe_ie<-hess_EM_dir_K1_oe(theta_u=theta.curr[t,], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess_re<-hess_EM_dir_K1_re(theta_u=theta.curr[t,], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess_oe_ie_re<-hess_EM_dir_K1_oe_re(theta_u=theta.curr[t,], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess<-matrix(c(hess_oe_ie,hess_oe_ie_re,hess_oe_ie_re,hess_re),2,2)
      #print(kappa(hess))
      theta.next[t,]<-theta.curr[t,]-(solve(hess)%*%gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    theta<-array(NA_real_,dim=c(T_grid,2,n_iter))
    theta[,,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,,iter_index]<-theta.update_K1(theta.curr=theta[,,iter_index-1], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids)
      #print(theta[,,iter_index])
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_dir_K1(theta_u = theta[t,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  ## Defining Model Selection functions for K=1 based on converged paramters.
  
  ## Defining a function to calculate the estimate of H_K
  H_K1<-function(theta,N,T_data){
    ## Assuming grid spacing is same as the data spacing.
    H_K_mat<-matrix(0,2,2)
    for (t in 1:T_data){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val_1<-exp(theta[t,1])
          exp_val_2<-exp(2*theta[t,2])
          H_K_mat[1,1]<-H_K_mat[1,1]+((2*exp_val_1*(1+exp_val_2))/((1+2*exp_val_1+exp_val_2)^2))
          H_K_mat[2,2]<-H_K_mat[2,2]+((4*exp_val_2*(1+2*exp_val_1))/((1+2*exp_val_1+exp_val_2)^2))
          H_K_mat[1,2]<-H_K_mat[1,2]-((4*exp_val_1*exp_val_2)/((1+2*exp_val_1+exp_val_2)^2))
        }
      }
    }
    H_K_mat[2,1]<-H_K_mat[1,2]
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K1<-function(theta,network,N,T_data){
    V_K_mat<-matrix(0,2,2)
    ## Assuming grid spacing is same as the data spacing.
    for (t in 1:T_data){
      u_vec<-rep(0,2)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val_1<-exp(theta[t,1])
          exp_val_2<-exp(2*theta[t,2])
          indicator_10<-(network[i,j,t]==1)&(network[j,i,t]==0)
          indicator_01<-(network[i,j,t]==0)&(network[j,i,t]==1)
          indicator_11<-(network[i,j,t]==1)&(network[j,i,t]==1)
          u_vec[1]<-u_vec[1]+(indicator_10+indicator_01-(2*exp_val_1/(1+2*exp_val_1+exp_val_2)))
          u_vec[2]<-u_vec[2]+(2*(indicator_11-(exp_val_2/(1+2*exp_val_1+exp_val_2))))
        }
      }
      V_K_mat<-V_K_mat+(u_vec%*%t(u_vec))
    }
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_data){
    ## Assuming grid spacing is same as the data spacing.
    cl_val<-0
    for(t in 1:T_data){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val_1<-exp(theta[t,1])
          exp_val_2<-exp(2*theta[t,2])
          indicator_10<-(network[i,j,t]==1)&(network[j,i,t]==0)
          indicator_01<-(network[i,j,t]==0)&(network[j,i,t]==1)
          indicator_11<-(network[i,j,t]==1)&(network[j,i,t]==1)
          cl_val<-cl_val+(indicator_10*theta[t,1]+indicator_01*theta[t,1]+2*indicator_11*theta[t,2]-(log(1+2*exp_val_1+exp_val_2)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC_K1<-function(theta,network,N,T_data,bandwidth){
    H_K_mat<-H_K1(theta = theta,N = N,T_data = T_data)
    V_K_mat<-V_K1(theta = theta,network = network, N = N,T_data = T_data)
    d_K<-sum(diag(solve(H_K_mat)%*%V_K_mat))
    cond_loglik_K1_val<-(cond_loglik_K1(theta = theta,network = network,N = N,T_data = T_data))
    t1<-(-2*cond_loglik_K1_val)
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(cond_loglik_K1_val,t2,BIC_val))
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE function
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the parameters 
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  #T_grid<-T_data
  
  ## Defining the actual data points and grid points 
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  ## Defining the number of clusters
  K<-nclust 
  #################################################
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
    for(i in 1:N){
      gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
      gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  if(K==1){
    start<-matrix(0,T_grid,2)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
  }else{
    start<-list()
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K)
    start[[3]]<-array(0,dim=c(T_grid,K,2))
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,,n_last]
    BIC_val<-BIC_K1(theta = param_converge,network = sim.net,N = N,T_data = T_data,bandwidth = h)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = rep(1,N),cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        RASE_theta_oe<-RASE_theta(theta_est = param_converge[,1],theta_true = theta_true[,1])
        RASE_theta_re<-RASE_theta(theta_est = param_converge[,2],theta_true = theta_true[,2])
        RASE_vec<-c(RASE_theta_oe,RASE_theta_re)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_vec)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    BIC_val<-BIC(theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data,bandwidth = h,cluster_ids=cluster_ids_est)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_oe_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est_oe<-param_converge[[3]][,K_permute_mat[k,],1]
          RASE_theta_oe_vec[k]<-RASE_theta(theta_est = theta_est_oe,theta_true = theta_true[,,1])
        }
        permute_true_id<-which.min(RASE_theta_oe_vec)
        permute_true<-K_permute_mat[permute_true_id,]
        theta_est_re<-param_converge[[3]][,permute_true,2]
        RASE_theta_oe<-RASE_theta_oe_vec[permute_true_id]
        RASE_theta_re<-RASE_theta(theta_est = theta_est_re,theta_true = theta_true[,,2])
        RASE_vec<-c(RASE_theta_oe,RASE_theta_re)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_vec)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper cross validation function to run for fixed number of clusters and a given bandwidth
wrapper_EM_dir_CV<-function(sim.net,nclust,bandwidth,fold_list,folds){
  ## This is the combined final updated and working code for NCV for EM directed case for all K
  #########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(EMdir)
  library(EMdirK1)
  library(Matrix)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids,test_node_set_len){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    quad_lin_coeff<-gamma_update_dir_CV(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
    
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of the theta with dimensions T_grid*K*2
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
    theta.next<-array(NA_real_,dim=c(T_grid,K,2))
    for (t in 1:T_grid){
      kernel_vec<-rep(0,T_data)
      for (temp_index in 1:T_data){
        kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      nonzero_ids<-which(kernel_vec!=0)
      nonzero_ids_len<-length(nonzero_ids)
      grad_vec_kern_oe<-matrix(0,T_data,K)
      grad_vec_kern_re<-matrix(0,T_data,K)
      if (nonzero_ids_len!=0){
        grad_oe<-grad_EM_dir_oe_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,test_node_set_len = test_node_set_len)
        grad_re<-grad_EM_dir_re_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,test_node_set_len = test_node_set_len)
        for (temp_index in 1:nonzero_ids_len){
          grad_vec_kern_oe[nonzero_ids[temp_index],]<-grad_oe[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
          grad_vec_kern_re[nonzero_ids[temp_index],]<-grad_re[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
        }
      }
      gradient_oe<-as.vector(colSums(grad_vec_kern_oe))
      gradient_re<-as.vector(colSums(grad_vec_kern_re))
      gradient<-c(gradient_oe,gradient_re)
      hess_oe<-hess_EM_dir_oe_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
      hess_re<-hess_EM_dir_re_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
      hess_oe_re<-hess_EM_dir_oe_re_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
      hess<-matrix(NA_real_,2*K,2*K)
      hess[1:K,1:K]<-hess_oe
      hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
      hess[(1:K),((K+1):(2*K))]<-hess_oe_re
      hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
      #theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)
      if(sum(hess)==0){
        theta.next[t,,]<-theta.curr[t,,]
      } else{theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)}
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the indices for actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(T_grid,K,2,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)&(iter_index<25)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta array
      theta[,,,iter_index]<-theta.update(theta.curr=theta[,,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len = test_node_set_len)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_dir_CV(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  ##########################################################################################################################################################################################################
  ## Defining a function to calculate the estimate of conditional log likelihood
    cond_loglik<-function(theta,network,N,K,T_data,cluster_ids,test_node_set_len){
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    cl_val<-0
    for(t in 1:T_data){
      for(i in 1:test_node_set_len){
        if(i!=test_node_set_len){
          for(j in (i+1):test_node_set_len){
            exp_val_1<-exp(theta[t,cluster_ids[i],1])
            exp_val_2<-exp(theta[t,cluster_ids[j],1])
            exp_val_3<-exp(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2])
            indicator_10<-(network[i,j,t]==1)&(network[i,j,t]==0)
            indicator_01<-(network[i,j,t]==0)&(network[i,j,t]==1)
            indicator_11<-(network[i,j,t]==1)&(network[i,j,t]==1)
            cl_val<-cl_val+((indicator_10*theta[t,cluster_ids[i],1])+(indicator_01*theta[t,cluster_ids[j],1])+(indicator_11*(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2]))-(log(1+exp_val_1+exp_val_2+exp_val_3)))
          }
        } else if(i==test_node_set_len){
          exp_val_1<-exp(theta[t,cluster_ids[i],1])
          exp_val_2<-exp(theta[t,cluster_ids[test_node_set_len],1])
          exp_val_3<-exp(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[test_node_set_len],2])
          cl_val<-cl_val+(-(log(1+exp_val_1+exp_val_2+exp_val_3)))
        }
      }
    }
    return(cl_val)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining the functions for K=1
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
    theta.next<-rep(NA_real_,T_grid)
    for (t in 1:T_grid){
      grad_vec<-grad_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data,test_node_set_len=test_node_set_len)
      grad_vec_kern<-rep(NA_real_,T_data)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-sum(grad_vec_kern)
      hess<-hess_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
      theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-matrix(NA_real_,N,n_iter)
    pi<-rep(NA_real_,n_iter)
    theta<-matrix(NA_real_,T_grid,n_iter)
    gamma[,1]<-start[[1]]
    pi[1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      #ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,iter_index]<-rep(1,N)
      ## Updating the pi vector
      pi[iter_index]<-1
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], pi=pi[iter_index],gamma=gamma[,iter_index], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len=test_node_set_len) 
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_K1(gamma=gamma[,iter_index], pi=pi[iter_index], theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      #print(iter_index)
      #print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  ##########################################################################################################################################################################################################
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids,test_node_set_len){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:test_node_set_len){
        if(i!=test_node_set_len){
          for(j in (i+1):test_node_set_len){
            exp_val<-exp(2*theta_data[t])
            cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
          }
        } else if(i==test_node_set_len){
          exp_val<-exp(2*theta_data[t])
          cl_val<-cl_val+((network[i,test_node_set_len,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  #################################################
  ## Defining the parameters 
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  #T_grid<-T_data
  
  ## Defining the actual data points and grid points 
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  K<-nclust ## Defining the number of clusters
  
  ## Initializing the cond_loglik_val vector to store the loss corresponding to each fold
  cond_loglik_val<-rep(NA_real_,folds)
  ## Loop over different folds for cross validation to choose bandwidth
  for (fold_index in 1:folds){
    test_node_set<-sort(fold_list[[fold_index]])
    train_list<-fold_list
    train_list[[fold_index]]<-NULL
    train_node_set<-sort(unlist(train_list))
    permute_set<-c(test_node_set,train_node_set)
    sim.net_permute<-array(NA_real_,dim = dim(sim.net))
    for (i in 1:T_data){
      sim.net_permute[,,i]<-as(as.integer(permute_set), "pMatrix")%*%sim.net[,,i]%*%t(as(as.integer(permute_set), "pMatrix"))
    }
    
    #################################################
    ## Setting the initial value of gamma
    if(K==1){
      gamma.start <- rep(1,N)
    }else{
      ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
      set.seed((2))
      MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
      gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
      ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
      for(i in 1:N){
        gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
        gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
      }
    }
    
    #################################################
    ## Defining the starting values of the iterator and running the main algorithm
    if(K==1){
      start<-matrix(0,T_grid,2)
      param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.00001,bandwidth = h, grid_ids = grid_ids)
    }else{
      start<-list()
      start[[1]]<-gamma.start
      start[[2]]<-rep(1/K,K)
      start[[3]]<-array(0,dim=c(T_grid,K,2))
      #debug(iterator)
      param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.001,bandwidth = h, grid_ids = grid_ids,test_node_set_len = length(test_node_set))
    }
    
    #################################################
    ## extracting the coverged parameter values and calculating BIC
    n_iter=1
    indicator_last<-0
    while(indicator_last==0){
      if(K==1){
        temp<-is.na(param[1,1,n_iter])
      }else{temp<-is.na(param[[1]][1,1,n_iter])}
      if(temp==TRUE){
        n_last<-n_iter-1
        indicator_last<-1
      }
      n_iter<-n_iter+1
    }
    param_converge<-list()
    if(K==1){
      param_converge[[1]]<-param[[1]][,n_last]
      param_converge[[2]]<-param[[2]][n_last]
      param_converge[[3]]<-param[[3]][,n_last]
      cond_loglik_val[fold_index]<-cond_loglik_K1(theta = param_converge[[3]],network = sim.net_permute,N = N,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = length(test_node_set))
    }else{
      param_converge[[1]]<-param[[1]][,,n_last]
      param_converge[[2]]<-param[[2]][,n_last]
      param_converge[[3]]<-param[[3]][,,,n_last]
      cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
          cluster_id<-which.max(param_converge[[1]][x,])
          return(cluster_id)
        }))
      cond_loglik_val[fold_index]<-cond_loglik(theta = param_converge[[3]],network = sim.net_permute,N = N,K = K,T_data = T_data,cluster_ids=cluster_ids,test_node_set_len = length(test_node_set))
    }
  }
  cond_loglik_sum<--sum(cond_loglik_val)
  return(cond_loglik_sum)
}

#########################################################################################################
## Defining a function band_wrapper for doing cross validation and choosing bandwidth for a given repeat
band_wrapper_EM_dir<-function(sim.net,bandwidths,K,folds){
  ###################################################
  ## Defining a wrapper cross validation function to run for fixed number of clusters and a given bandwidth
  wrapper_EM_dir_CV<-function(sim.net,nclust,bandwidth,fold_list,folds){
    ## This is the combined final updated and working code for NCV for EM directed case for all K
    #########################################################################################################
    ## Loading the required packages
    require(lda)
    library(quadprog)
    library(combinat)
    library(EMdir)
    library(EMdirK1)
    library(Matrix)
    
    ######################################################################################################
    ## Defining the epanechnikov kernel function
    epan<-function(y){
      return(0.75*(1-y^2)*(abs(y)<=1))
    }
    
    #################################################
    ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
    ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
    ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
    gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids,test_node_set_len){
      gamma.next<-matrix(NA_real_,N,K)
      constraint_matrix<-matrix(NA_real_,2+K,K)
      constraint_matrix[1,]<-rep(1,K)
      constraint_matrix[2,]<-rep(-1,K)
      constraint_matrix[3:(K+2),]<-diag(K)
      constraint_vector<-c(1,-1,rep(0,K))
      quad_lin_coeff<-gamma_update_dir_CV(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
      
      for (i in 1:N){
        gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
      }
      
      ## normalizing gamma_i. deal with case outside (0,1) later
      gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
        x_norm<-x/sum(x)
        return(x_norm)
      }))
      return(gamma.next)
    }
    
    #################################################
    ## Defining a function to update K*1 vector of pi
    pi.update<-function(gamma.curr,N,K){
      pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
      ## normalization of pi
      pi.next<-pi.next/sum(pi.next)
      return(pi.next)
    }
    
    #################################################
    ## Defining the update function of the theta with dimensions T_grid*K*2
    theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
      theta.next<-array(NA_real_,dim=c(T_grid,K,2))
      for (t in 1:T_grid){
        kernel_vec<-rep(0,T_data)
        for (temp_index in 1:T_data){
          kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
        }
        nonzero_ids<-which(kernel_vec!=0)
        nonzero_ids_len<-length(nonzero_ids)
        grad_vec_kern_oe<-matrix(0,T_data,K)
        grad_vec_kern_re<-matrix(0,T_data,K)
        if (nonzero_ids_len!=0){
          grad_oe<-grad_EM_dir_oe_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,test_node_set_len = test_node_set_len)
          grad_re<-grad_EM_dir_re_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,test_node_set_len = test_node_set_len)
          for (temp_index in 1:nonzero_ids_len){
            grad_vec_kern_oe[nonzero_ids[temp_index],]<-grad_oe[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
            grad_vec_kern_re[nonzero_ids[temp_index],]<-grad_re[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
          }
        }
        gradient_oe<-as.vector(colSums(grad_vec_kern_oe))
        gradient_re<-as.vector(colSums(grad_vec_kern_re))
        gradient<-c(gradient_oe,gradient_re)
        hess_oe<-hess_EM_dir_oe_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
        hess_re<-hess_EM_dir_re_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
        hess_oe_re<-hess_EM_dir_oe_re_CV(theta_u=as.matrix(theta.curr[t,,]), gamma=gamma, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
        hess<-matrix(NA_real_,2*K,2*K)
        hess[1:K,1:K]<-hess_oe
        hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
        hess[(1:K),((K+1):(2*K))]<-hess_oe_re
        hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
        #theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)
        if(sum(hess)==0){
          theta.next[t,,]<-theta.curr[t,,]
        } else{theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)}
      }
      return(theta.next)
    }
    
    #################################################
    ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
    iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
      N<-dim(network)[1]
      T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
      ## adjacency matrix
      T_grid<-length(grid_ids)  ## Number of grid points
      
      ## Defining the indices for actual data points
      data_ids<-1:T_data
      
      ## initializing the arrays for parameters
      gamma<-array(NA_real_,dim=c(N,K,n_iter))
      pi<-matrix(NA_real_,K,n_iter)
      theta<-array(NA_real_,dim=c(T_grid,K,2,n_iter))
      gamma[,,1]<-start[[1]]
      pi[,1]<-start[[2]]
      theta[,,,1]<-start[[3]]
      ## iterations
      ## Defining the iteration index
      iter_index<-2
      ## Initializing the error
      error<-Inf
      ## Initializing the current ELBO values over the whole grid
      ELBO_grid.curr<-rep(10^10,T_grid)
      while((error>thres)&(iter_index<25)){
        ## Starting the stopwatch to calculate the per iteration time
        ptm<-proc.time()
        ## Updating the N*K gamma matrix i.e. variational variational parameters
        gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = test_node_set_len)
        ## Updating the pi vector
        pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
        ## Updating the theta array
        theta[,,,iter_index]<-theta.update(theta.curr=theta[,,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len = test_node_set_len)
        ELBO_grid.prev<-ELBO_grid.curr
        ELBO_grid.curr<-rep(NA_real_,T_grid)
        for(t in 1:T_grid){
          ELBO_grid.curr[t]<-ELBO_conv_EM_dir_CV(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len = test_node_set_len)
        }
        error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
        print(error)
        print(iter_index)
        print(proc.time()-ptm)
        iter_index<-iter_index+1
      }
      return(list(gamma,pi,theta))
    }
    
    ##########################################################################################################################################################################################################
    ## Defining a function to calculate the estimate of conditional log likelihood
    cond_loglik<-function(theta,network,N,K,T_data,cluster_ids,test_node_set_len){
      ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
      cl_val<-0
      for(t in 1:T_data){
        for(i in 1:test_node_set_len){
          if(i!=test_node_set_len){
            for(j in (i+1):test_node_set_len){
              exp_val_1<-exp(theta[t,cluster_ids[i],1])
              exp_val_2<-exp(theta[t,cluster_ids[j],1])
              exp_val_3<-exp(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2])
              indicator_10<-(network[i,j,t]==1)&(network[i,j,t]==0)
              indicator_01<-(network[i,j,t]==0)&(network[i,j,t]==1)
              indicator_11<-(network[i,j,t]==1)&(network[i,j,t]==1)
              cl_val<-cl_val+((indicator_10*theta[t,cluster_ids[i],1])+(indicator_01*theta[t,cluster_ids[j],1])+(indicator_11*(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2]))-(log(1+exp_val_1+exp_val_2+exp_val_3)))
            }
          } else if(i==test_node_set_len){
            exp_val_1<-exp(theta[t,cluster_ids[i],1])
            exp_val_2<-exp(theta[t,cluster_ids[test_node_set_len],1])
            exp_val_3<-exp(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[test_node_set_len],2])
            cl_val<-cl_val+(-(log(1+exp_val_1+exp_val_2+exp_val_3)))
          }
        }
      }
      return(cl_val)
    }
    
    #########################################################################################################
    #########################################################################################################
    ## Defining the functions for K=1
    #################################################
    ## Defining the update function of full matrix theta with dimensions T_grid*K
    theta.update_K1<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,test_node_set_len){
      theta.next<-rep(NA_real_,T_grid)
      for (t in 1:T_grid){
        grad_vec<-grad_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data,test_node_set_len=test_node_set_len)
        grad_vec_kern<-rep(NA_real_,T_data)
        for (temp_index in 1:T_data){
          grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
        }
        gradient<-sum(grad_vec_kern)
        hess<-hess_ELBO_K1_CV(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
        theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
      }
      return(theta.next)
    }
    
    #################################################
    ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
    iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids,test_node_set_len){
      N<-dim(network)[1]
      T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
      ## adjacency matrix
      T_grid<-length(grid_ids)  ## Number of grid points
      
      ## Defining the actual data points
      data_ids<-1:T_data
      
      ## initializing the arrays for parameters
      gamma<-matrix(NA_real_,N,n_iter)
      pi<-rep(NA_real_,n_iter)
      theta<-matrix(NA_real_,T_grid,n_iter)
      gamma[,1]<-start[[1]]
      pi[1]<-start[[2]]
      theta[,1]<-start[[3]]
      ## iterations
      ## Defining the iteration index
      iter_index<-2
      ## Initializing the error
      error<-Inf
      ## Initializing the current ELBO values over the whole grid
      ELBO_grid.curr<-rep(10^10,T_grid)
      while(error>thres){
        ## Starting the stopwatch to calculate the per iteration time
        #ptm<-proc.time()
        ## Updating the N*K gamma matrix i.e. variational variational parameters
        gamma[,iter_index]<-rep(1,N)
        ## Updating the pi vector
        pi[iter_index]<-1
        ## Updating the theta matrix
        theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], pi=pi[iter_index],gamma=gamma[,iter_index], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,test_node_set_len=test_node_set_len) 
        ELBO_grid.prev<-ELBO_grid.curr
        ELBO_grid.curr<-rep(NA_real_,T_grid)
        for(t in 1:T_grid){
          ELBO_grid.curr[t]<-ELBO_conv_K1(gamma=gamma[,iter_index], pi=pi[iter_index], theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),test_node_set_len=test_node_set_len)
        }
        error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
        #print(iter_index)
        #print(proc.time()-ptm)
        iter_index<-iter_index+1
      }
      return(list(gamma,pi,theta))
    }
    
    ##########################################################################################################################################################################################################
    ## Defining a function to calculate the estimate of conditional log likelihood
    cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids,test_node_set_len){
      ## Assuming grid spacing is half the data spacing.
      theta_data<-theta
      cl_val<-0
      for(t in 1:T_grid){
        for(i in 1:test_node_set_len){
          if(i!=test_node_set_len){
            for(j in (i+1):test_node_set_len){
              exp_val<-exp(2*theta_data[t])
              cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
            }
          } else if(i==test_node_set_len){
            exp_val<-exp(2*theta_data[t])
            cl_val<-cl_val+((network[i,test_node_set_len,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
          }
        }
      }
      return(cl_val)
    }
    
    #################################################
    ## Defining the parameters 
    N<-dim(sim.net)[1] ## Number of nodes from the network
    ## Total time points in the network for which we have data in the form of adjacency matrix
    T_data<-dim(sim.net)[3]
    
    #T_grid<-T_data
    
    ## Defining the actual data points and grid points 
    data_ids<-(1:T_data)
    #grid_ids<-1:T_data
    grid_ids<-seq(1,T_data,by = 1)
    
    ## Defining the total number of grid points
    T_grid<-length(grid_ids)
    
    # Defining the epanechnikov kernel bandwidth: h
    h<-bandwidth
    
    K<-nclust ## Defining the number of clusters
    
    ## Initializing the cond_loglik_val vector to store the loss corresponding to each fold
    cond_loglik_val<-rep(NA_real_,folds)
    ## Loop over different folds for cross validation to choose bandwidth
    for (fold_index in 1:folds){
      test_node_set<-sort(fold_list[[fold_index]])
      train_list<-fold_list
      train_list[[fold_index]]<-NULL
      train_node_set<-sort(unlist(train_list))
      permute_set<-c(test_node_set,train_node_set)
      sim.net_permute<-array(NA_real_,dim = dim(sim.net))
      for (i in 1:T_data){
        sim.net_permute[,,i]<-as(as.integer(permute_set), "pMatrix")%*%sim.net[,,i]%*%t(as(as.integer(permute_set), "pMatrix"))
      }
      
      #################################################
      ## Setting the initial value of gamma
      if(K==1){
        gamma.start <- rep(1,N)
      }else{
        ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
        set.seed((2))
        MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
        gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
        ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
        for(i in 1:N){
          gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
          gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
        }
      }
      
      #################################################
      ## Defining the starting values of the iterator and running the main algorithm
      if(K==1){
        start<-matrix(0,T_grid,2)
        param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.00001,bandwidth = h, grid_ids = grid_ids)
      }else{
        start<-list()
        start[[1]]<-gamma.start
        start[[2]]<-rep(1/K,K)
        start[[3]]<-array(0,dim=c(T_grid,K,2))
        #debug(iterator)
        param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.001,bandwidth = h, grid_ids = grid_ids,test_node_set_len = length(test_node_set))
      }
      
      #################################################
      ## extracting the coverged parameter values and calculating BIC
      n_iter=1
      indicator_last<-0
      while(indicator_last==0){
        if(K==1){
          temp<-is.na(param[1,1,n_iter])
        }else{temp<-is.na(param[[1]][1,1,n_iter])}
        if(temp==TRUE){
          n_last<-n_iter-1
          indicator_last<-1
        }
        n_iter<-n_iter+1
      }
      param_converge<-list()
      if(K==1){
        param_converge[[1]]<-param[[1]][,n_last]
        param_converge[[2]]<-param[[2]][n_last]
        param_converge[[3]]<-param[[3]][,n_last]
        cond_loglik_val[fold_index]<-cond_loglik_K1(theta = param_converge[[3]],network = sim.net_permute,N = N,T_grid = T_grid,grid_ids=grid_ids,test_node_set_len = length(test_node_set))
      }else{
        param_converge[[1]]<-param[[1]][,,n_last]
        param_converge[[2]]<-param[[2]][,n_last]
        param_converge[[3]]<-param[[3]][,,,n_last]
        cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
          cluster_id<-which.max(param_converge[[1]][x,])
          return(cluster_id)
        }))
        cond_loglik_val[fold_index]<-cond_loglik(theta = param_converge[[3]],network = sim.net_permute,N = N,K = K,T_data = T_data,cluster_ids=cluster_ids,test_node_set_len = length(test_node_set))
      }
    }
    cond_loglik_sum<--sum(cond_loglik_val)
    return(cond_loglik_sum)
  }
  
  ###################################################
  ## Defining a function to get N_Vfolds used as a input to wrapper_CV type of functions while doing NCV
  get_fold_list<-function(sim.net,folds){
    ## Dividing the nodes set into V folds equal subsets where V = 3 is most commonly used
    ## randomly sampling the nodes set
    N<-dim(sim.net)[1] ## Number of nodes from the network
    N_set_rand<-sample(N)
    N_Vfolds<-list()
    step_size<-floor(N/folds)
    for (i in 1:folds){
      N_Vfolds[[i]]<-N_set_rand[(1+(i-1)*step_size):(i*step_size)]
    }
    if((N%%folds)!=0){
      for (i in 1:(N%%folds)){
        N_Vfolds[[i]]<-c(N_Vfolds[[i]],N_set_rand[i+(step_size*folds)])
      }
    }
    return(N_Vfolds)
  }
  
  ###################################################
  ## Initializing the conditional log-likelihood vector for different bandwidths
  cond_loglik_vec<-rep(NA_real_,length(bandwidths))
  for (band_index in 1:length(bandwidths)){
    cond_loglik_vec[band_index]<- wrapper_EM_dir_CV(sim.net = sim.net,nclust = K,bandwidth = bandwidths[band_index],fold_list = get_fold_list(sim.net = sim.net,folds = folds),folds = folds)
  }
  bandwidth_chosen<-bandwidths[which.min(cond_loglik_vec)]
  return(bandwidth_chosen)
}

#########################################################################################################
## Defining a wrapper function for SEM undirected case to run for different number of clusters and bandwidths
wrapper_SEM_undir<-function(sim.net,nclust,bandwidth){
  ## This is the combined final updated and working code for SEM undirected case for all K
  #######################################################################################################
  ######################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(SEMundir)
  #library(SEMundirK1)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N.data*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a N.grid*K matrix for current theta estimates which is a subset of full 
  ## matrix theta. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,cluster_ids){
    theta.next<-matrix(NA_real_,T_grid,K)
    for (t in 1:T_grid){
      kernel_vec<-rep(0,T_data)
      for (temp_index in 1:T_data){
        kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      nonzero_ids<-which(kernel_vec!=0)
      nonzero_ids_len<-length(nonzero_ids)
      grad_vec_kern<-matrix(0,T_data,K)
      if (nonzero_ids_len!=0){
        grad_vec<-grad_SEM_undir(theta_u=as.vector(theta.curr[t,]), network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,cluster_ids=cluster_ids)
        for (temp_index in 1:nonzero_ids_len){
          grad_vec_kern[nonzero_ids[temp_index],]<-grad_vec[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
        }
      }
      gradient<-colSums(grad_vec_kern)
      hess<-hess_SEM_undir(theta_u=as.vector(theta.curr[t,]), N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),cluster_ids=cluster_ids)
      if(sum(hess)==0){
        theta.next[t,]<-as.vector(theta.curr[t,])
      } else{theta.next[t,]<-as.vector(theta.curr[t,])-as.vector(solve(hess)%*%gradient)}
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the indices for actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(T_grid,K,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)&(iter_index<25)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      #debug(gamma.update.wrapper)
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Stochastic step
      ## Sampling the cluster memberships for each node
      cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        y<-as.vector(rmultinom(n = 1,size = 1,prob = gamma[x,,iter_index]))
        cluster_id<-which(y==1)
        return(cluster_id)
      }))
      ## Updating the theta matrix
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,cluster_ids=cluster_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_SEM_undir(gamma=gamma[,,iter_index],pi=pi[,iter_index], theta_u = theta[t,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),cluster_ids = cluster_ids)
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Defining Model Selection functions based on converged paramters for K>=2.
  
  ## Defining a function to calculate the estimate of H_K
  H_K<-function(gamma,theta,N,K,T_grid){
    H_K_mat<-matrix(NA_real_,K,K)
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    ## Filling the upper triangular matrix of H_K since it is symmetric
    for (k in 1:K){
      for (l in k:K){
        if (k==l){
          t1<-0
          for (t  in 1:T_grid){
            for(i in 1:N){
              if(i!=N){
                for(j in (i+1):N){
                  indicator_1<-gamma[i,k]>(1/K)
                  indicator_2<-gamma[j,k]>(1/K)
                  if(indicator_1|indicator_2){
                    indicator<-indicator_1+indicator_2
                    exp_val<-exp(2*theta_data[t,k])
                    t1<-t1+((exp_val/((1+exp_val)^2))*(indicator^2))
                  }
                }
              } else if(i==N){
                indicator_1<-gamma[i,k]>(1/K)
                indicator_2<-gamma[N,k]>(1/K)
                if(indicator_1|indicator_2){
                  indicator<-indicator_1+indicator_2
                  exp_val<-exp(2*theta_data[t,k])
                  t1<-t1+((exp_val/((1+exp_val)^2))*(indicator^2))
                }
              }
            }
          }
          H_K_mat[k,k]<-t1
        }else if (k!=l){
          t1<-0
          for (t  in 1:T_grid){
            for(i in 1:N){
              if(i!=N){
                for(j in (i+1):N){
                  if(((gamma[i,k]>(1/K))&(gamma[j,l]>(1/K)))|((gamma[i,l]>(1/K))&(gamma[j,k]>(1/K)))){
                    exp_val<-exp(theta_data[t,k]+theta_data[t,l])
                    t1<-t1+(exp_val/((1+exp_val)^2))
                  }
                }
              } else if(i==N){
                if(((gamma[i,k]>(1/K))&(gamma[N,l]>(1/K)))|((gamma[i,l]>(1/K))&(gamma[N,k]>(1/K)))){
                  exp_val<-exp(theta_data[t,k]+theta_data[t,l])
                  t1<-t1+(exp_val/((1+exp_val)^2))
                }
              }
            }
          }
          H_K_mat[k,l]<-t1
          H_K_mat[l,k]<-t1
        }
      }
    }
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K<-function(gamma,theta,network,N,K,T_grid,grid_ids){
    V_K_mat<-matrix(0,K,K)
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    for (t in 1:T_grid){
      u_vec<-rep(0,K)
      for (k in 1:K){
        for(i in 1:N){
          if(i!=N){
            for(j in (i+1):N){
              indicator_i<-(gamma[i,k]>(1/K))
              indicator_j<-(gamma[j,k]>(1/K))
              if(indicator_i|indicator_j){
                if(indicator_i&indicator_j){
                  exp_val<-exp(2*theta_data[t,k])
                  u_vec[k]<-u_vec[k]+((network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))*(indicator_i+indicator_j))
                }else if(indicator_i&(!indicator_j)){
                  clusterid_j<-which.max(gamma[j,])
                  exp_val<-exp(theta_data[t,k]+theta_data[t,clusterid_j])
                  u_vec[k]<-u_vec[k]+(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))
                }else if(indicator_j&(!indicator_i)){
                  clusterid_i<-which.max(gamma[i,])
                  exp_val<-exp(theta_data[t,k]+theta_data[t,clusterid_i])
                  u_vec[k]<-u_vec[k]+(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))
                }
              }
            }
          } else if(i==N){
            indicator_i<-(gamma[i,k]>(1/K))
            indicator_j<-(gamma[N,k]>(1/K))
            if(indicator_i|indicator_j){
              if(indicator_i&indicator_j){
                exp_val<-exp(2*theta_data[t,k])
                u_vec[k]<-u_vec[k]+((network[i,N,grid_ids[t]]-(exp_val/(1+exp_val)))*(indicator_i+indicator_j))
              }else if(indicator_i&(!indicator_j)){
                clusterid_j<-which.max(gamma[N,])
                exp_val<-exp(theta_data[t,k]+theta_data[t,clusterid_j])
                u_vec[k]<-u_vec[k]+(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))
              }else if(indicator_j&(!indicator_i)){
                clusterid_i<-which.max(gamma[i,])
                exp_val<-exp(theta_data[t,k]+theta_data[t,clusterid_i])
                u_vec[k]<-u_vec[k]+(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val)))
              }
            }
          }
        }
      }
      V_K_mat<-V_K_mat+u_vec%*%t(u_vec)
    }
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(gamma,theta,network,N,K,T_grid,grid_ids){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-which.max(gamma[i,])
          cluster_id_j<-which.max(gamma[j,])
          exp_val<-exp(theta_data[t,cluster_id_i]+theta_data[t,cluster_id_j])
          cl_val<-cl_val+((network[i,j,grid_ids[t]]*(theta_data[t,cluster_id_i]+theta_data[t,cluster_id_j]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC<-function(gamma,theta,network,N,K,T_data,T_grid,bandwidth,grid_ids){
    H_K_mat<-H_K(gamma = gamma,theta = theta,N = N,K = K,T_grid = T_grid)
    V_K_mat<-V_K(gamma = gamma,theta = theta,network = network, N = N,K = K,T_grid = T_grid,grid_ids=grid_ids)
    d_K<-sum(diag(solve(H_K_mat+(10^(-10)*diag(K)))%*%(V_K_mat)))
    #d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    t1<-(-2*(cond_loglik(gamma = gamma,theta = theta,network = network,N = N,K = K,T_grid = T_grid,grid_ids=grid_ids)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  #######################################################################################################
  #######################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-rep(NA_real_,T_grid)
    for (t in 1:T_grid){
      grad_vec<-grad_ELBO_K1(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data)
      grad_vec_kern<-rep(NA_real_,T_data)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-sum(grad_vec_kern)
      hess<-hess_ELBO_K1(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-matrix(NA_real_,N,n_iter)
    pi<-rep(NA_real_,n_iter)
    theta<-matrix(NA_real_,T_grid,n_iter)
    gamma[,1]<-start[[1]]
    pi[1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      #ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,iter_index]<-rep(1,N)
      ## Updating the pi vector
      pi[iter_index]<-1
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], pi=pi[iter_index],gamma=gamma[,iter_index], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids) 
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_K1(gamma=gamma[,iter_index], pi=pi[iter_index], theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      #print(iter_index)
      #print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Defining Model Selection functions for K=1 based on converged paramters.
  
  ## Defining a function to calculate the estimate of H_K
  H_K1<-function(theta,N,T_grid){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    H_K_val<-0
    for (t in 1:T_grid){
      for(i in 1:N){
        if(i!=N){
          for(j in (i+1):N){
            exp_val<-exp(2*theta_data[t])
            H_K_val<-H_K_val+(4*(exp_val/((1+exp_val)^2)))
          }
        } else if(i==N){
          exp_val<-exp(2*theta_data[t])
          H_K_val<-H_K_val+(4*(exp_val/((1+exp_val)^2)))
        }
      }
    }
    return(H_K_val)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K1<-function(theta,network,N,T_grid,grid_ids){
    V_K_val<-0
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    for (t in 1:T_grid){
      u_val<-0
      for(i in 1:N){
        if(i!=N){
          for(j in (i+1):N){
            exp_val<-exp(2*theta_data[t])
            u_val<-u_val+(2*(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val))))
          }
        } else if(i==N){
          exp_val<-exp(2*theta_data[t])
          u_val<-u_val+(2*(network[i,N,grid_ids[t]]-(exp_val/(1+exp_val))))
        }
      }
      V_K_val<-V_K_val+(u_val^2)
    }
    return(V_K_val)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC_K1<-function(gamma,theta,network,N,K,T_data,T_grid,bandwidth,grid_ids){
    H_K_val<-H_K1(theta = theta,N = N,T_grid = T_grid)
    V_K_val<-V_K1(theta = theta,network = network, N = N,T_grid = T_grid,grid_ids=grid_ids)
    d_K<-V_K_val/H_K_val
    t1<-(-2*(cond_loglik_K1(theta = theta,network = network,N = N,T_grid = T_grid,grid_ids=grid_ids)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  #######################################################################################################
  ## Defining the parameters 
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  #T_grid<-T_data
  
  ## Defining the actual data points and grid points 
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  K<-nclust ## Defining the number of clusters
  
  #################################################
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
    for(i in 1:N){
      gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
      gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  #undebug(iterator_K1)
  if(K==1){
    start[[1]]<-gamma.start
    start[[2]]<-1
    start[[3]]<-rep(0,T_grid)
    #debug(iterator_K1)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.001,bandwidth = h, grid_ids = grid_ids)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K)#pi.update(gamma.curr = gamma.start,N=N,K=K)
    start[[3]]<-matrix(0,T_grid,K)
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.001,bandwidth = h, grid_ids = grid_ids)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[[1]][1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge[[1]]<-param[[1]][,n_last]
    param_converge[[2]]<-param[[2]][n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    BIC_val<-BIC_K1(theta = param_converge[[3]],network = sim.net,N = N,T_data = T_data,T_grid = T_grid,bandwidth = h,grid_ids=grid_ids)
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,n_last]
    BIC_val<-BIC(gamma = param_converge[[1]],theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data,T_grid = T_grid,bandwidth = h,grid_ids=grid_ids)
  }
  output_list<-list(param_converge,BIC_val)
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for SEM undirected case to run for different number of clusters and bandwidths
wrapper_SEM_dir<-function(sim.net,nclust,bandwidth,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM directed case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(SEMdir)
  library(EMdirK1)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_grid,grid_ids){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    quad_lin_coeff<-gamma_update_dir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K, T_grid=T_grid,grid_ids=grid_ids)
    
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of the theta with dimensions T_grid*K*2
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids,cluster_ids){
    theta.next<-array(NA_real_,dim=c(T_grid,K,2))
    for (t in 1:T_grid){
      kernel_vec<-rep(0,T_data)
      for (temp_index in 1:T_data){
        kernel_vec[temp_index]<-(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      nonzero_ids<-which(kernel_vec!=0)
      nonzero_ids_len<-length(nonzero_ids)
      grad_vec_kern_oe<-matrix(0,T_data,K)
      grad_vec_kern_re<-matrix(0,T_data,K)
      if (nonzero_ids_len!=0){
        grad_oe<-grad_SEM_dir_oe(theta_u=as.matrix(theta.curr[t,,]), network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,cluster_ids=cluster_ids)
        grad_re<-grad_SEM_dir_re(theta_u=as.matrix(theta.curr[t,,]), network=network, N=N, K=K, nonzero_ids=nonzero_ids,nonzero_ids_len=nonzero_ids_len,cluster_ids=cluster_ids)
        for (temp_index in 1:nonzero_ids_len){
          grad_vec_kern_oe[nonzero_ids[temp_index],]<-grad_oe[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
          grad_vec_kern_re[nonzero_ids[temp_index],]<-grad_re[temp_index,]*kernel_vec[nonzero_ids[temp_index]]
        }
      }
      gradient_oe<-as.vector(colSums(grad_vec_kern_oe))
      gradient_re<-as.vector(colSums(grad_vec_kern_re))
      gradient<-c(gradient_oe,gradient_re)
      hess_oe<-hess_SEM_dir_oe(theta_u=as.matrix(theta.curr[t,,]), N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),cluster_ids=cluster_ids)
      hess_re<-hess_SEM_dir_re(theta_u=as.matrix(theta.curr[t,,]), N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),cluster_ids=cluster_ids)
      hess_oe_re<-hess_SEM_dir_oe_re(theta_u=as.matrix(theta.curr[t,,]), N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1),cluster_ids=cluster_ids)
      hess<-matrix(NA_real_,2*K,2*K)
      hess[1:K,1:K]<-hess_oe
      hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
      hess[(1:K),((K+1):(2*K))]<-hess_oe_re
      hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
      #theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)
      if(sum(hess)==0){
        theta.next[t,,]<-theta.curr[t,,]
      } else{theta.next[t,,]<-matrix(c(theta.curr[t,,1],theta.curr[t,,2])-as.vector(solve(hess)%*%gradient),K,2)}
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the indices for actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(T_grid,K,2,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)&(iter_index<30)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,,iter_index-1],network=network,N=N,K=K,T_grid = T_grid,grid_ids=grid_ids)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Stochastic step
      ## Sampling the cluster memberships for each node
      cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        y<-as.vector(rmultinom(n = 1,size = 1,prob = gamma[x,,iter_index]))
        cluster_id<-which(y==1)
        return(cluster_id)
      }))
      ## Updating the theta array
      theta[,,,iter_index]<-theta.update(theta.curr=theta[,,,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data=T_data,T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids,cluster_ids = cluster_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_SEM_dir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta_u = theta[t,,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1), cluster_ids=cluster_ids)
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  ######################################################################################################
  #################################################
  ## Model Selection functions for K>=2
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(theta,network,N,K,T_data,cluster_ids){
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    cl_val<-0
    for(t in 1:T_data){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val_1<-exp(theta[t,cluster_ids[i],1])
          exp_val_2<-exp(theta[t,cluster_ids[j],1])
          exp_val_3<-exp(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2])
          indicator_10<-(network[i,j,t]==1)&(network[j,i,t]==0)
          indicator_01<-(network[i,j,t]==0)&(network[j,i,t]==1)
          indicator_11<-(network[i,j,t]==1)&(network[j,i,t]==1)
          cl_val<-cl_val+((indicator_10*theta[t,cluster_ids[i],1])+(indicator_01*theta[t,cluster_ids[j],1])+(indicator_11*(theta[t,cluster_ids[i],2]+theta[t,cluster_ids[j],2]))-(log(1+exp_val_1+exp_val_2+exp_val_3)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of H_K
  H_K<-function(theta,N,K,T_data,cluster_ids){
    H_K_mat<-matrix(NA_real_,2*K,2*K)
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    H_K_oe<-HK_dir_oe(theta=theta, N=N, K=K, T_data=T_data, cluster_ids=cluster_ids)
    H_K_re<-HK_dir_re(theta=theta, N=N, K=K, T_data=T_data, cluster_ids=cluster_ids)
    H_K_oe_re<-HK_dir_oe_re(theta=theta, N=N, K=K, T_data=T_data, cluster_ids=cluster_ids)
    H_K_mat[1:K,1:K]<-H_K_oe
    H_K_mat[((K+1):(2*K)),((K+1):(2*K))]<-H_K_re
    H_K_mat[(1:K),((K+1):(2*K))]<-H_K_oe_re
    H_K_mat[((K+1):(2*K)),(1:K)]<-t(H_K_oe_re)
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K<-function(theta,network,N,K,T_data,cluster_ids){
    V_K_mat<-matrix(0,2*K,2*K)
    ## Assuming grid spacing is same as the data spacing i.e. grid points coincide with data points
    u_mat_oe<-VK_dir_oe(theta=theta, network=network, N=N, K=K, T_data=T_data,cluster_ids=cluster_ids)
    u_mat_re<-VK_dir_re(theta=theta, network=network, N=N, K=K, T_data=T_data,cluster_ids=cluster_ids)
    for (t in 1:T_data){
      u_vec<-c(as.vector(u_mat_oe[t,]),as.vector(u_mat_re[t,]))
      V_K_mat<-V_K_mat+u_vec%*%t(u_vec)
    }
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC<-function(theta,network,N,K,T_data,bandwidth,cluster_ids){
    H_K_mat<-H_K(theta = theta,N = N,K = K,T_data = T_data,cluster_ids = cluster_ids)
    V_K_mat<-V_K(theta = theta,network = network, N = N,K = K,T_data = T_data,cluster_ids=cluster_ids)
    d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    cond_loglik_val<-cond_loglik(theta = theta,network = network,N = N,K = K,T_data = T_data,cluster_ids=cluster_ids)
    t1<-(-2*cond_loglik_val)
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(cond_loglik_val,t2,BIC_val))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of theta_oe/ie with dimensions T_grid
  theta.update_K1<-function(theta.curr,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-matrix(NA_real_,T_grid,2)
    for (t in 1:T_grid){
      grad_oe_ie<-grad_EM_dir_K1_oe(theta_u=theta.curr[t,], network=network, N=N, K=K, T_data=T_data)
      grad_re<-grad_EM_dir_K1_re(theta_u=theta.curr[t,], network=network, N=N, K=K, T_data=T_data)
      grad_vec_kern<-matrix(NA_real_,T_data,2)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index,1]<-grad_oe_ie[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
        grad_vec_kern[temp_index,2]<-grad_re[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-as.vector(colSums(grad_vec_kern))
      hess_oe_ie<-hess_EM_dir_K1_oe(theta_u=theta.curr[t,], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess_re<-hess_EM_dir_K1_re(theta_u=theta.curr[t,], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess_oe_ie_re<-hess_EM_dir_K1_oe_re(theta_u=theta.curr[t,], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      hess<-matrix(c(hess_oe_ie,hess_oe_ie_re,hess_oe_ie_re,hess_re),2,2)
      #print(kappa(hess))
      theta.next[t,]<-theta.curr[t,]-(solve(hess)%*%gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    theta<-array(NA_real_,dim=c(T_grid,2,n_iter))
    theta[,,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while((error>thres)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,,iter_index]<-theta.update_K1(theta.curr=theta[,,iter_index-1], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids)
      #print(theta[,,iter_index])
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_dir_K1(theta_u = theta[t,,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  ## Defining Model Selection functions for K=1 based on converged paramters.
  
  ## Defining a function to calculate the estimate of H_K
  H_K1<-function(theta,N,T_data){
    ## Assuming grid spacing is same as the data spacing.
    H_K_mat<-matrix(0,2,2)
    for (t in 1:T_data){
      for(i in 1:N){
        if(i!=N){
          for(j in (i+1):N){
            exp_val_1<-exp(theta[t,1])
            exp_val_2<-exp(2*theta[t,2])
            H_K_mat[1,1]<-H_K_mat[1,1]+((2*exp_val_1*(1+exp_val_2))/((1+2*exp_val_1+exp_val_2)^2))
            H_K_mat[2,2]<-H_K_mat[2,2]+((4*exp_val_2*(1+2*exp_val_1))/((1+2*exp_val_1+exp_val_2)^2))
            H_K_mat[1,2]<-H_K_mat[1,2]-((4*exp_val_1*exp_val_2)/((1+2*exp_val_1+exp_val_2)^2))
          }
        } else if(i==N){
          exp_val_1<-exp(theta[t,1])
          exp_val_2<-exp(2*theta[t,2])
          H_K_mat[1,1]<-H_K_mat[1,1]+((2*exp_val_1*(1+exp_val_2))/((1+2*exp_val_1+exp_val_2)^2))
          H_K_mat[2,2]<-H_K_mat[2,2]+((4*exp_val_2*(1+2*exp_val_1))/((1+2*exp_val_1+exp_val_2)^2))
          H_K_mat[1,2]<-H_K_mat[1,2]-((4*exp_val_1*exp_val_2)/((1+2*exp_val_1+exp_val_2)^2))
        }
      }
    }
    H_K_mat[2,1]<-H_K_mat[1,2]
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K1<-function(theta,network,N,T_data){
    V_K_mat<-matrix(0,2,2)
    ## Assuming grid spacing is same as the data spacing.
    for (t in 1:T_data){
      u_vec<-rep(0,2)
      for(i in 1:N){
        if(i!=N){
          for(j in (i+1):N){
            exp_val_1<-exp(theta[t,1])
            exp_val_2<-exp(2*theta[t,2])
            indicator_10<-(network[i,j,t]==1)&(network[j,i,t]==0)
            indicator_01<-(network[i,j,t]==0)&(network[j,i,t]==1)
            indicator_11<-(network[i,j,t]==1)&(network[j,i,t]==1)
            u_vec[1]<-u_vec[1]+(indicator_10+indicator_01-(2*exp_val_1/(1+2*exp_val_1+exp_val_2)))
            u_vec[2]<-u_vec[2]+(2*(indicator_11-(exp_val_2/(1+2*exp_val_1+exp_val_2))))
          }
        } else if(i==N){
          exp_val_1<-exp(theta[t,1])
          exp_val_2<-exp(2*theta[t,2])
          u_vec[1]<-u_vec[1]+(-(2*exp_val_1/(1+2*exp_val_1+exp_val_2)))
          u_vec[2]<-u_vec[2]+(2*(-(exp_val_2/(1+2*exp_val_1+exp_val_2))))
        }
      }
      V_K_mat<-V_K_mat+(u_vec%*%t(u_vec))
    }
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_data){
    ## Assuming grid spacing is same as the data spacing.
    cl_val<-0
    for(t in 1:T_data){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val_1<-exp(theta[t,1])
          exp_val_2<-exp(2*theta[t,2])
          indicator_10<-(network[i,j,t]==1)&(network[j,i,t]==0)
          indicator_01<-(network[i,j,t]==0)&(network[j,i,t]==1)
          indicator_11<-(network[i,j,t]==1)&(network[j,i,t]==1)
          cl_val<-cl_val+(indicator_10*theta[t,1]+indicator_01*theta[t,1]+2*indicator_11*theta[t,2]-(log(1+2*exp_val_1+exp_val_2)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC_K1<-function(theta,network,N,T_data,bandwidth){
    H_K_mat<-H_K1(theta = theta,N = N,T_data = T_data)
    V_K_mat<-V_K1(theta = theta,network = network, N = N,T_data = T_data)
    d_K<-sum(diag(solve(H_K_mat)%*%V_K_mat))
    cond_loglik_K1_val<-(cond_loglik_K1(theta = theta,network = network,N = N,T_data = T_data))
    t1<-(-2*cond_loglik_K1_val)
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(cond_loglik_K1_val,t2,BIC_val))
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE function
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the parameters 
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  #T_grid<-T_data
  
  ## Defining the actual data points and grid points 
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  ## Defining the number of clusters
  K<-nclust 
  #################################################
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 1000)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
    for(i in 1:N){
      gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
      gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  if(K==1){
    start<-matrix(0,T_grid,2)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
  }else{
    start<-list()
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K)
    start[[3]]<-array(0,dim=c(T_grid,K,2))
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,,n_last]
    BIC_val<-BIC_K1(theta = param_converge,network = sim.net,N = N,T_data = T_data,bandwidth = h)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = rep(1,N),cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        RASE_theta_oe<-RASE_theta(theta_est = param_converge[,1],theta_true = theta_true[,1])
        RASE_theta_re<-RASE_theta(theta_est = param_converge[,2],theta_true = theta_true[,2])
        RASE_vec<-c(RASE_theta_oe,RASE_theta_re)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_vec)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    BIC_val<-BIC(theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data,bandwidth = h,cluster_ids=cluster_ids_est)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_oe_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est_oe<-param_converge[[3]][,K_permute_mat[k,],1]
          RASE_theta_oe_vec[k]<-RASE_theta(theta_est = theta_est_oe,theta_true = theta_true[,,1])
        }
        permute_true_id<-which.min(RASE_theta_oe_vec)
        permute_true<-K_permute_mat[permute_true_id,]
        theta_est_re<-param_converge[[3]][,permute_true,2]
        RASE_theta_oe<-RASE_theta_oe_vec[permute_true_id]
        RASE_theta_re<-RASE_theta(theta_est = theta_est_re,theta_true = theta_true[,,2])
        RASE_vec<-c(RASE_theta_oe,RASE_theta_re)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_vec)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }
  return(output_list)
}

#########################################################################################################
wrapper_sc_EM_undir<-function(sim.net,nclust,thres,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(scEMundir)
  library(EMundirK1)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_sc_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_sc_EM_undir(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K)
    hess<-hess_sc_EM_undir(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,iter_index-1],network=network,N=N,K=K)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_sc_EM_undir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta = theta[,iter_index], network=network, N=N, K=K)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Defining Model Selection functions based on converged paramters for K>=2.
  
  ## Defining a function to calculate the estimate of H_K
  H_K<-function(gamma,theta,N,K){
    H_K_mat<-matrix(NA_real_,K,K)
    ## Filling the upper triangular matrix of H_K since it is symmetric
    for (k in 1:K){
      for (l in k:K){
        if (k==l){
          t1<-0
          for(i in 1:(N-1)){
            for(j in (i+1):N){
              indicator_1<-gamma[i,k]>(1/K)
              indicator_2<-gamma[j,k]>(1/K)
              if(indicator_1|indicator_2){
                indicator<-indicator_1+indicator_2
                exp_val<-exp(2*theta[k])
                t1<-t1+((exp_val/((1+exp_val)^2))*(indicator^2))
              }
            }
          }
          H_K_mat[k,k]<-t1
        }else if (k!=l){
          t1<-0
          for(i in 1:(N-1)){
            for(j in (i+1):N){
              if(((gamma[i,k]>(1/K))&(gamma[j,l]>(1/K)))|((gamma[i,l]>(1/K))&(gamma[j,k]>(1/K)))){
                exp_val<-exp(theta[k]+theta[l])
                t1<-t1+(exp_val/((1+exp_val)^2))
              }
            }
          }
          H_K_mat[k,l]<-t1
          H_K_mat[l,k]<-t1
        }
      }
    }
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K<-function(gamma,theta,network,N,K){
    V_K_mat<-matrix(0,K,K)
    u_vec<-rep(0,K)
    for (k in 1:K){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          indicator_i<-(gamma[i,k]>(1/K))
          indicator_j<-(gamma[j,k]>(1/K))
          if(indicator_i|indicator_j){
            if(indicator_i&indicator_j){
              exp_val<-exp(2*theta[k])
              u_vec[k]<-u_vec[k]+((network[i,j]-(exp_val/(1+exp_val)))*(indicator_i+indicator_j))
            }else if(indicator_i&(!indicator_j)){
              clusterid_j<-which.max(gamma[j,])
              exp_val<-exp(theta[k]+theta[clusterid_j])
              u_vec[k]<-u_vec[k]+(network[i,j]-(exp_val/(1+exp_val)))
            }else if(indicator_j&(!indicator_i)){
              clusterid_i<-which.max(gamma[i,])
              exp_val<-exp(theta[k]+theta[clusterid_i])
              u_vec[k]<-u_vec[k]+(network[i,j]-(exp_val/(1+exp_val)))
            }
          }
        }
      }
    }
    V_K_mat<-V_K_mat+u_vec%*%t(u_vec)
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(gamma,theta,network,N,K){
    cl_val<-0
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        cluster_id_i<-which.max(gamma[i,])
        cluster_id_j<-which.max(gamma[j,])
        exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
        cl_val<-cl_val+((network[i,j]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC<-function(gamma,theta,network,N,K){
    H_K_mat<-H_K(gamma = gamma,theta = theta,N = N,K = K)
    V_K_mat<-V_K(gamma = gamma,theta = theta,network = network, N = N,K = K)
    d_K<-sum(diag(solve(H_K_mat+(10^(-10)*diag(K)))%*%(V_K_mat)))
    #d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    t1<-(-2*(cond_loglik(gamma = gamma,theta = theta,network = network,N = N,K = K)))
    t2<-d_K*log((N*(N-1))/2)*2.1153*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-rep(NA_real_,T_grid)
    for (t in 1:T_grid){
      grad_vec<-grad_EM_undir_K1(theta_u=theta.curr[t], network=network, N=N, K=K, T_data=T_data)
      grad_vec_kern<-rep(NA_real_,T_data)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-sum(grad_vec_kern)
      hess<-hess_EM_undir_K1(theta_u=theta.curr[t], N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    theta<-matrix(NA_real_,T_grid,n_iter)
    theta[,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_EM_undir_K1(theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  ## Defining Model Selection functions for K=1 based on converged paramters.
  
  ## Defining a function to calculate the estimate of H_K
  H_K1<-function(theta,N,T_grid){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    H_K_val<-0
    for (t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          H_K_val<-H_K_val+(4*(exp_val/((1+exp_val)^2)))
        }
      }
    }
    return(H_K_val)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K1<-function(theta,network,N,T_grid,grid_ids){
    V_K_val<-0
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    for (t in 1:T_grid){
      u_val<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          u_val<-u_val+(2*(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val))))
        }
      }
      V_K_val<-V_K_val+(u_val^2)
    }
    return(V_K_val)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC_K1<-function(theta,network,N,K,T_data,T_grid,bandwidth,grid_ids){
    H_K_val<-H_K1(theta = theta,N = N,T_grid = T_grid)
    V_K_val<-V_K1(theta = theta,network = network, N = N,T_grid = T_grid,grid_ids=grid_ids)
    d_K<-V_K_val/H_K_val
    t1<-(-2*(cond_loglik_K1(theta = theta,network = network,N = N,T_grid = T_grid,grid_ids=grid_ids)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE function
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  
  K<-nclust ## Defining the number of clusters
  
  #################################################
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net,alpha = 1/2,num.iterations = 100)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
    for(i in 1:N){
      gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
      gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-0
    #debug(iterator_K1)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres,bandwidth = h, grid_ids = grid_ids)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K)#pi.update(gamma.curr = gamma.start,N=N,K=K)
    start[[3]]<-rep(0,K)
    #debug(iterator)
    #ptm<-proc.time()
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
    #print(proc.time()-ptm)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,n_last]
    BIC_val<-BIC_K1(theta = param_converge,network = sim.net,N = N,T_data = T_data,T_grid = T_grid,bandwidth = h,grid_ids=grid_ids)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = rep(1,N),cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_theta)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    BIC_val<-BIC(gamma = param_converge[[1]], theta = param_converge[[3]],network = sim.net,N = N,K = K)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][,K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id<-which.min(RASE_theta_vec)
        RASE_theta<-RASE_theta_vec[permute_true_id]
        output_list<-list(param_converge,BIC_val,RI_val,RASE_theta)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for SEM undirected case to run for different number of clusters and bandwidths
wrapper_sc_SEM_undir<-function(sim.net,nclust,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for SEM undirected case for all K
  #######################################################################################################
  ######################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(scSEMundir)
  
  ######################################################################################################
  ## Defining the epanechnikov kernel function
  epan<-function(y){
    return(0.75*(1-y^2)*(abs(y)<=1))
  }
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N.data*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a N.grid*K matrix for current theta estimates which is a subset of full 
  ## matrix theta. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    quad_lin_coeff<-gamma_update_sc_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K)
    #print(quad_lin_coeff)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full vector theta with dimensions K*1
  theta.update<-function(theta.curr,network,N,K,cluster_ids){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_sc_SEM_undir(theta=theta.curr, network=network, N=N, K=K, cluster_ids=cluster_ids)
    hess<-hess_sc_SEM_undir(theta=theta.curr, N=N, K=K,cluster_ids=cluster_ids)
    #print(hess)
    theta.next<-theta.curr-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<20)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      #debug(gamma.update.wrapper)
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,iter_index-1],network=network,N=N,K=K)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Stochastic step
      ## Sampling the cluster memberships for each node
      cluster_ids<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        y<-as.vector(rmultinom(n = 1,size = 1,prob = gamma[x,,iter_index]))
        cluster_id<-which(y==1)
        return(cluster_id)
      }))
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], network=network, N=N, K=K, cluster_ids=cluster_ids)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_sc_SEM_undir(gamma=gamma[,,iter_index],pi=pi[,iter_index], theta = theta[,iter_index], network=network, N=N, K=K, cluster_ids = cluster_ids)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Defining Model Selection functions based on converged paramters for K>=2.
  ## Defining a function to calculate the estimate of H_K
  H_K<-function(gamma,theta,N,K){
    H_K_mat<-matrix(NA_real_,K,K)
    ## Filling the upper triangular matrix of H_K since it is symmetric
    for (k in 1:K){
      for (l in k:K){
        if (k==l){
          t1<-0
          for(i in 1:(N-1)){
            for(j in (i+1):N){
              indicator_1<-gamma[i,k]>(1/K)
              indicator_2<-gamma[j,k]>(1/K)
              if(indicator_1|indicator_2){
                indicator<-indicator_1+indicator_2
                exp_val<-exp(2*theta[k])
                t1<-t1+((exp_val/((1+exp_val)^2))*(indicator^2))
              }
            }
          }
          H_K_mat[k,k]<-t1
        }else if (k!=l){
          t1<-0
          for(i in 1:(N-1)){
            for(j in (i+1):N){
              if(((gamma[i,k]>(1/K))&(gamma[j,l]>(1/K)))|((gamma[i,l]>(1/K))&(gamma[j,k]>(1/K)))){
                exp_val<-exp(theta[k]+theta[l])
                t1<-t1+(exp_val/((1+exp_val)^2))
              }
            }
          }
          H_K_mat[k,l]<-t1
          H_K_mat[l,k]<-t1
        }
      }
    }
    return(H_K_mat)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K<-function(gamma,theta,network,N,K){
    u_vec<-rep(0,K)
    for (k in 1:K){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          indicator_i<-(gamma[i,k]>(1/K))
          indicator_j<-(gamma[j,k]>(1/K))
          if(indicator_i|indicator_j){
            if(indicator_i&indicator_j){
              exp_val<-exp(2*theta[k])
              u_vec[k]<-u_vec[k]+((network[i,j]-(exp_val/(1+exp_val)))*(indicator_i+indicator_j))
            }else if(indicator_i&(!indicator_j)){
              clusterid_j<-which.max(gamma[j,])
              exp_val<-exp(theta[k]+theta[clusterid_j])
              u_vec[k]<-u_vec[k]+(network[i,j]-(exp_val/(1+exp_val)))
            }else if(indicator_j&(!indicator_i)){
              clusterid_i<-which.max(gamma[i,])
              exp_val<-exp(theta[k]+theta[clusterid_i])
              u_vec[k]<-u_vec[k]+(network[i,j]-(exp_val/(1+exp_val)))
            }
          }
        }
      }
    }
    V_K_mat<-u_vec%*%t(u_vec)
    return(V_K_mat)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik<-function(gamma,theta,network,N,K){
    cl_val<-0
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        cluster_id_i<-which.max(gamma[i,])
        cluster_id_j<-which.max(gamma[j,])
        exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
        cl_val<-cl_val+((network[i,j]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC<-function(gamma,theta,network,N,K){
    H_K_mat<-H_K(gamma = gamma,theta = theta,N = N,K = K)
    V_K_mat<-V_K(gamma = gamma,theta = theta,network = network, N = N,K = K)
    d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    #d_K<-sum(diag(solve(H_K_mat)%*%(V_K_mat)))
    t1<-(-2*(cond_loglik(gamma = gamma,theta = theta,network = network,N = N,K = K)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  #######################################################################################################
  #######################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,pi,gamma,network,N,K,T_data,T_grid,bandwidth,data_ids,grid_ids){
    theta.next<-rep(NA_real_,T_grid)
    for (t in 1:T_grid){
      grad_vec<-grad_ELBO_K1(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data)
      grad_vec_kern<-rep(NA_real_,T_data)
      for (temp_index in 1:T_data){
        grad_vec_kern[temp_index]<-grad_vec[temp_index]*(1/bandwidth)*epan((data_ids[temp_index]-grid_ids[t])/bandwidth)
      }
      gradient<-sum(grad_vec_kern)
      hess<-hess_ELBO_K1(theta_u=theta.curr[t], gamma=gamma, network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      theta.next[t]<-theta.curr[t]-((1/hess)*gradient)
    }
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres,bandwidth,grid_ids){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of 
    ## adjacency matrix
    T_grid<-length(grid_ids)  ## Number of grid points
    
    ## Defining the actual data points
    data_ids<-1:T_data
    
    ## initializing the arrays for parameters
    gamma<-matrix(NA_real_,N,n_iter)
    pi<-rep(NA_real_,n_iter)
    theta<-matrix(NA_real_,T_grid,n_iter)
    gamma[,1]<-start[[1]]
    pi[1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-rep(10^10,T_grid)
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      #ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,iter_index]<-rep(1,N)
      ## Updating the pi vector
      pi[iter_index]<-1
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], pi=pi[iter_index],gamma=gamma[,iter_index], network=network, N=N, K=K,T_data=T_data, T_grid = T_grid, bandwidth=bandwidth, data_ids=data_ids,grid_ids = grid_ids) 
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-rep(NA_real_,T_grid)
      for(t in 1:T_grid){
        ELBO_grid.curr[t]<-ELBO_conv_K1(gamma=gamma[,iter_index], pi=pi[iter_index], theta_u = theta[t,iter_index], network=network, N=N, K=K, T_data=T_data, bandwidth=bandwidth, data_ids=data_ids, grid_ids=grid_ids, grid_id_index=(t-1))
      }
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      #print(iter_index)
      #print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  #################################################
  ## Defining Model Selection functions for K=1 based on converged paramters.
  
  ## Defining a function to calculate the estimate of H_K
  H_K1<-function(theta,N,T_grid){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    H_K_val<-0
    for (t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          H_K_val<-H_K_val+(4*(exp_val/((1+exp_val)^2)))
        }
      }
    }
    return(H_K_val)
  }
  
  ## Defining a function to calculate the estimate of V_K
  V_K1<-function(theta,network,N,T_grid,grid_ids){
    V_K_val<-0
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    for (t in 1:T_grid){
      u_val<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          u_val<-u_val+(2*(network[i,j,grid_ids[t]]-(exp_val/(1+exp_val))))
        }
      }
      V_K_val<-V_K_val+(u_val^2)
    }
    return(V_K_val)
  }
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  cond_loglik_K1<-function(theta,network,N,T_grid,grid_ids){
    ## Assuming grid spacing is half the data spacing.
    theta_data<-theta
    cl_val<-0
    for(t in 1:T_grid){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          exp_val<-exp(2*theta_data[t])
          cl_val<-cl_val+((network[i,j,grid_ids[t]]*(2*theta_data[t]))-(log(1+exp_val)))
        }
      }
    }
    return(cl_val)
  }
  
  ## Defining a function to calculate the estimate of BIC for a given K
  BIC_K1<-function(theta,network,N,K,T_data,T_grid,bandwidth,grid_ids){
    H_K_val<-H_K1(theta = theta,N = N,T_grid = T_grid)
    V_K_val<-V_K1(theta = theta,network = network, N = N,T_grid = T_grid,grid_ids=grid_ids)
    d_K<-V_K_val/H_K_val
    t1<-(-2*(cond_loglik_K1(theta = theta,network = network,N = N,T_grid = T_grid,grid_ids=grid_ids)))
    t2<-d_K*log(T_data*((N*(N-1))/2))*2.1153*(1/bandwidth)*T_data*0.45
    BIC_val<-t1+t2
    return(c(t1,t2,BIC_val))
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE function
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  ## Defining the actual data points and grid points
  data_ids<-(1:T_data)
  #grid_ids<-1:T_data
  grid_ids<-seq(1,T_data,by = 1)
  
  ## Defining the total number of grid points
  T_grid<-length(grid_ids)
  
  # Defining the epanechnikov kernel bandwidth: h
  h<-bandwidth
  
  K<-nclust ## Defining the number of clusters
  
  #################################################
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 10000)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
    for(i in 1:N){
      gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
      gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-0
    #debug(iterator_K1)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001,bandwidth = h, grid_ids = grid_ids)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K)#pi.update(gamma.curr = gamma.start,N=N,K=K)
    start[[3]]<-rep(0,K)
    #debug(iterator)
    #ptm<-proc.time()
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=0.0001)
    #print(proc.time()-ptm)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,n_last]
    BIC_val<-BIC_K1(theta = param_converge,network = sim.net,N = N,T_data = T_data,T_grid = T_grid,bandwidth = h,grid_ids=grid_ids)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = rep(1,N),cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,BIC_val,RI_val,RASE_theta)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    BIC_val<-BIC(gamma = param_converge[[1]], theta = param_converge[[3]],network = sim.net,N = N,K = K)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][,K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id<-which.min(RASE_theta_vec)
        RASE_theta<-RASE_theta_vec[permute_true_id]
        output_list<-list(param_converge,BIC_val,RI_val,RASE_theta)
      }else{
        output_list<-list(param_converge,BIC_val,RI_val)
      }
    }else{output_list<-list(param_converge,BIC_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a function to simulate the data in the form of a sequence of adjacency matrices from a undirected or directed network. Inputs: 
## 1) Number of nodes(N),
## 2) For Dynamic network, dyn=1, else 0
## 3) For Directed network, dir=1, else 0 
## 4) theta: For directed network: theta is array T_data*K*2 or theta is matrix T_data*2 for K=1; For undirected network: theta is matrix T_data*K or theta is vector T_data*1 for K=1 and 
## 5) pi(mixing proportions) (pi is 1 as default for K=1)
simulate_network<-function(N,dyn,dir,theta,pi=1){
  if (dyn==1){
    if(is.vector(theta)!=1){
      T_data<-dim(theta)[1]
    }else{T_data<-length(theta)}
    network<-array(NA_integer_,dim = c(N,N,T_data))
    if (length(pi)>1){
      clust_id_sampling<-rmultinom(n = N, size = 1, prob = pi)
      clust_ids<-apply(X = clust_id_sampling,MARGIN = 2,FUN = function(x){
        y<-which(x==1)
        return(y)
      })
      if(dir==1){
        for (t in 1:T_data){
          for (i in 1:(N-1)){
            for (j in (i+1):N){
              ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
              exp_val_1<-exp(theta[t,clust_ids[i],1])
              exp_val_2<-exp(theta[t,clust_ids[j],1])
              exp_val_3<-exp(theta[t,clust_ids[i],2]+theta[t,clust_ids[j],2])
              probab_00<-1/(1+exp_val_1+exp_val_2+exp_val_3)
              probab_10<-exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)
              probab_01<-exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)
              probab_11<-exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)
              probab<-c(probab_00,probab_10,probab_01,probab_11)
              #print(probab)
              edge_sample<-as.vector(rmultinom(n = 1, size = 1, prob = probab))
              edge_id<-which(edge_sample==1)
              if(edge_id==1){
                network[i,j,t]<-0
                network[j,i,t]<-0
              }else if(edge_id==2){
                network[i,j,t]<-1
                network[j,i,t]<-0
              }else if(edge_id==3){
                network[i,j,t]<-0
                network[j,i,t]<-1
              }else if(edge_id==4){
                network[i,j,t]<-1
                network[j,i,t]<-1
              }
            }
          }
          diag(network[,,t])<-0
        }
      }else if(dir==0){
        for (t in 1:T_data){
          for (i in 1:(N-1)){
            for (j in (i+1):N){
              ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
              exp_val<-exp(theta[t,clust_ids[i]]+theta[t,clust_ids[j]])
              probab<-exp_val/(1+exp_val)
              #print(1-probab,probab)
              edge_sample<-rbinom(n = 1, size = 1, prob = probab)
              if(edge_sample==0){
                network[i,j,t]<-0
                network[j,i,t]<-0
              }else if(edge_sample==1){
                network[i,j,t]<-1
                network[j,i,t]<-1
              }
            }
          }
          diag(network[,,t])<-0
        }
      }
      
    }else if(pi==1){
      clust_ids<-rep(1,N)
      if(dir==1){
        for (t in 1:T_data){
          for (i in 1:(N-1)){
            for (j in (i+1):N){
              ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
              exp_val_1<-exp(theta[t,1])
              exp_val_2<-exp(theta[t,1])
              exp_val_3<-exp(2*theta[t,2])
              probab_00<-1/(1+exp_val_1+exp_val_2+exp_val_3)
              probab_10<-exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)
              probab_01<-exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)
              probab_11<-exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)
              probab<-c(probab_00,probab_10,probab_01,probab_11)
              #print(probab)
              edge_sample<-as.vector(rmultinom(n = 1, size = 1, prob = probab))
              edge_id<-which(edge_sample==1)
              if(edge_id==1){
                network[i,j,t]<-0
                network[j,i,t]<-0
              }else if(edge_id==2){
                network[i,j,t]<-1
                network[j,i,t]<-0
              }else if(edge_id==3){
                network[i,j,t]<-0
                network[j,i,t]<-1
              }else if(edge_id==4){
                network[i,j,t]<-1
                network[j,i,t]<-1
              }
            }
          }
          diag(network[,,t])<-0
        }
      }else if(dir==0){
        for (t in 1:T_data){
          for (i in 1:(N-1)){
            for (j in (i+1):N){
              ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
              exp_val<-exp(2*theta[t])
              probab<-exp_val/(1+exp_val)
              #print(1-probab,probab)
              edge_sample<-rbinom(n = 1, size = 1, prob = probab)
              if(edge_sample==0){
                network[i,j,t]<-0
                network[j,i,t]<-0
              }else if(edge_sample==1){
                network[i,j,t]<-1
                network[j,i,t]<-1
              }
            }
          }
          diag(network[,,t])<-0
        }
      }
    }
  }else if(dyn==0){
    network<-matrix(NA_integer_,N,N)
    if (length(pi)>1){
      clust_id_sampling<-rmultinom(n = N, size = 1, prob = pi)
      clust_ids<-apply(X = clust_id_sampling,MARGIN = 2,FUN = function(x){
        y<-which(x==1)
        return(y)
      })
      if(dir==1){
        for (i in 1:(N-1)){
          for (j in (i+1):N){
            ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
            exp_val_1<-exp(theta[clust_ids[i],1])
            exp_val_2<-exp(theta[clust_ids[j],1])
            exp_val_3<-exp(theta[clust_ids[i],2]+theta[clust_ids[j],2])
            probab_00<-1/(1+exp_val_1+exp_val_2+exp_val_3)
            probab_10<-exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)
            probab_01<-exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)
            probab_11<-exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)
            probab<-c(probab_00,probab_10,probab_01,probab_11)
            #print(probab)
            edge_sample<-as.vector(rmultinom(n = 1, size = 1, prob = probab))
            edge_id<-which(edge_sample==1)
            if(edge_id==1){
              network[i,j]<-0
              network[j,i]<-0
            }else if(edge_id==2){
              network[i,j]<-1
              network[j,i]<-0
            }else if(edge_id==3){
              network[i,j]<-0
              network[j,i]<-1
            }else if(edge_id==4){
              network[i,j]<-1
              network[j,i]<-1
            }
          }
        }
        diag(network)<-0
      }else if(dir==0){
        for (i in 1:(N-1)){
          for (j in (i+1):N){
            ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
            exp_val<-exp(theta[clust_ids[i]]+theta[clust_ids[j]])
            probab<-exp_val/(1+exp_val)
            #print(1-probab,probab)
            edge_sample<-rbinom(n = 1, size = 1, prob = probab)
            if(edge_sample==0){
              network[i,j]<-0
              network[j,i]<-0
            }else if(edge_sample==1){
              network[i,j]<-1
              network[j,i]<-1
            }
          }
        }
        diag(network)<-0
      }
      
    }else if(pi==1){
      clust_ids<-rep(1,N)
      if(dir==1){
        for (i in 1:(N-1)){
          for (j in (i+1):N){
            ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
            exp_val_1<-exp(theta[1])
            exp_val_2<-exp(theta[1])
            exp_val_3<-exp(2*theta[2])
            probab_00<-1/(1+exp_val_1+exp_val_2+exp_val_3)
            probab_10<-exp_val_1/(1+exp_val_1+exp_val_2+exp_val_3)
            probab_01<-exp_val_2/(1+exp_val_1+exp_val_2+exp_val_3)
            probab_11<-exp_val_3/(1+exp_val_1+exp_val_2+exp_val_3)
            probab<-c(probab_00,probab_10,probab_01,probab_11)
            #print(probab)
            edge_sample<-as.vector(rmultinom(n = 1, size = 1, prob = probab))
            edge_id<-which(edge_sample==1)
            if(edge_id==1){
              network[i,j]<-0
              network[j,i]<-0
            }else if(edge_id==2){
              network[i,j]<-1
              network[j,i]<-0
            }else if(edge_id==3){
              network[i,j]<-0
              network[j,i]<-1
            }else if(edge_id==4){
              network[i,j]<-1
              network[j,i]<-1
            }
          }
        }
        diag(network)<-0
      }else if(dir==0){
        for (i in 1:(N-1)){
          for (j in (i+1):N){
            ## Defining the probabilities for D_i,j = (0,0), (1,0), (0,1), (1,1) given the cluster ids for node i and node j. Note D_ij=(1,0) implies Y(i,j)=1 and Y(j,i)=0 and vice versa.
            exp_val<-exp(2*theta)
            probab<-exp_val/(1+exp_val)
            #print(1-probab,probab)
            edge_sample<-rbinom(n = 1, size = 1, prob = probab)
            if(edge_sample==0){
              network[i,j]<-0
              network[j,i]<-0
            }else if(edge_sample==1){
              network[i,j]<-1
              network[j,i]<-1
            }
          }
        }
        diag(network)<-0
      }
    }
  }
  
  return(list(network,clust_ids))
}

#########################################################################################################
## Defining a function to get N_Vfolds used as a input to wrapper_CV type of functions while doing NCV
get_fold_list<-function(sim.net,folds){
  ## Dividing the nodes set into V folds equal subsets where V = 3 is most commonly used
  ## randomly sampling the nodes set
  N<-dim(sim.net)[1] ## Number of nodes from the network
  N_set_rand<-sample(N)
  N_Vfolds<-list()
  step_size<-floor(N/folds)
  for (i in 1:folds){
    N_Vfolds[[i]]<-N_set_rand[(1+(i-1)*step_size):(i*step_size)]
  }
  if((N%%folds)!=0){
    for (i in 1:(N%%folds)){
      N_Vfolds[[i]]<-c(N_Vfolds[[i]],N_set_rand[i+(step_size*folds)])
    }
  }
  return(N_Vfolds)
}

#########################################################################################################
## Creating simulation module for simulating the network for directed case K>=2 and testing the EM_undir code
sim_module_EM_undir<-function(){
  ## Defining theta matrix T*K undirected case
  T_data<-51 ## Defining number of data points
  K<-2 ## Defining the number of clusters
  z <- seq(0,1,by=0.02)
  theta_true<-matrix(NA_real_,T_data,K) 
  for(i in 1:length(z)){
    theta_true[i,1]<-(1/2)*cos(3*pi*z[i])
    theta_true[i,2]<--(1/2)*cos(3*pi*z[i])-1.75
    #theta_true[i,2]<-(1/2)*exp(z[i])+0.5
  }
  # 
  df_theta_true<-data.frame(1:51,theta_true)
  names(df_theta_true)<-c("time","K1","K2")
  df_theta_true<-data.frame(melt(df_theta_true,id.vars = "time"),id="True")
  p<-ggplot(data = df_theta_true)
  p_sim<-p+geom_point(mapping = aes(x = time,y = value,color=variable))+geom_line(mapping = aes(x = time,y = value,color=variable))+labs(x="time",y="theta")
  p_sim
  
  save(theta_true,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/theta_true/undir_K2_s1_theta_true.RData")
  #########################################################################################################
  
  # Simulating the network
  for (sim_index in 1:100){
    net_list<-simulate_network(N = 100,dyn = 1,dir = 0,theta = theta_true,pi = c(0.5,0.5))
    save(net_list,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/data/undir/K2/s1/n100/undir_K2_s1_n100_",sim_index,".RData"))
  }
  
  ## simulating the network
  net<-simulate_network(N = 100,dyn = 1,dir = 0,theta = theta_true,pi = c(0.5,0.5))
  
  #########################################################################################################
  ## Testing the EM undirected network code on the simulated dataset
  #undebug(wrapper_EM_undir)
  net_list<-wrapper_EM_undir(sim.net = net[[1]],nclust=3,bandwidth = 4,sim_indicator = 1,theta_true = theta_true,K_true = 2,cluster_ids_true = net[[2]])
  param<-net_list[[1]]
  
  ## Plotting for K=2
  df_theta_est<-data.frame(1:51,param[[3]])
  names(df_theta_est)<-c("time","K1","K2")
  df_theta_est<-data.frame(melt(df_theta_est,id.vars = "time"),id="Estimated")
  df_theta_final<-rbind(df_theta_est,df_theta_true)
  p<-ggplot(data = df_theta_final)
  p1<-p+geom_point(mapping = aes(x = time,y = value,linetype=id,color=variable),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=id,color=variable))+labs(x="Time",y="Parameter",linetype="",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
  p1
  return(p1)
}

#########################################################################################################
## Creating simulation module for simulating the network for directed case K>=2 and testing the EM_undir code
sim_module_EM_undir_K1<-function(){
  ## Defining theta matrix T*K undirected case
  T_data<-51 ## Defining number of data points
  K<-1 ## Defining the number of clusters
  z <- seq(0,1,by=0.02)
  theta_true<-matrix(NA_real_,T_data,K) 
  for(i in 1:length(z)){
    theta_true[i,1]<-cos(3*pi*z[i])
  }
  
  df_theta<-data.frame(1:51,theta_true)
  names(df_theta)<-c("time","K1")
  df_theta<-data.frame(melt(df_theta,id.vars = "time"),id="true")
  p<-ggplot(data = df_theta)
  p_sim<-p+geom_point(mapping = aes(x = time,y = value,color=variable))+geom_line(mapping = aes(x = time,y = value,color=variable))+labs(x="time",y="theta")
  p_sim
  
  save(theta_true,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/theta_true/undir_K1_s1_theta_true.RData")
  #########################################################################################################
  
  # Simulating the network
  for (sim_index in 1:100){
    net_list<-simulate_network(N = 100,dyn = 1,dir = 0,theta = theta_true,pi = 1)
    save(net_list,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/data/undir/K1/s1/n100/undir_K1_s1_n100_",sim_index,".RData"))
  }
  
  ## simulating the network
  net<-simulate_network(N = 100,dyn = 1,dir = 0,theta = as.vector(theta_true),pi = 1)
  
  #########################################################################################################
  ## Testing the EM undirected network code on the simulated dataset
  #undebug(wrapper_EM_undir)
  net_list<-wrapper_EM_undir(sim.net = net[[1]],nclust=1,bandwidth = 4,sim_indicator = 1,theta_true = theta_true,K_true = 1,cluster_ids_true = net[[2]])
  
  #BIC_results<-matrix(NA_real_,4,3)
  RASE<-rep(NA_real_,100)
  theta<-matrix(NA_real_,100,51)
  for (sim_index in 1:100){
    net<-simulate_network(N = 100,dyn = 1,dir = 0,theta = as.vector(theta_true),pi = 1)
    net_list<-wrapper_EM_undir(sim.net = net[[1]],nclust=1,bandwidth = 2.5,sim_indicator = 1,theta_true = theta_true,K_true = 1,cluster_ids_true = net[[2]])
    theta[sim_index,]<-net_list[[1]]
    RASE[sim_index]<-net_list[[4]]
    #BIC_results[k-1,]<-net_list[[2]]
  }
  RASE[order(RASE)][51]
  which(order(RASE)==51)
  theta_est<-theta[41,]
  
  ## Plotting for K=1
  df_theta<-data.frame("Time"=1:51,"Estimated"=theta_est,"True"=theta_true)
  df_theta_final<-melt(df_theta,id.vars = "Time")
  p<-ggplot(data = df_theta_final)
  p+geom_point(mapping = aes(x = Time,y = value,linetype=variable),size=0.2)+geom_line(mapping = aes(x = Time,y = value,linetype=variable))+labs(x="Time",y="Edge parameter",linetype="",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
  
  
  param<-net_list[[1]]
  df_theta_est<-data.frame(1:51,param[[3]])
  names(df_theta_est)<-c("time","K1","K2")
  df_theta_est<-data.frame(melt(df_theta_est,id.vars = "time"),id="Estimated")
  df_theta_final<-rbind(df_theta,df_theta_est)
  p<-ggplot(data = df_theta_final)
  p1<-p+geom_point(mapping = aes(x = time,y = value,linetype=variable,color=id),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=variable,color=id))+labs(x="time",y="theta")
  return(p1)
}

#########################################################################################################
## Creating simulation module for simulating the network for directed case K>=2 and testing the EM_undir code
sim_module_sc_EM_undir<-function(){
  ## Defining theta vector K*1 static undirected case
  ## simulating the network
  K<-11
  theta_true<-seq(from = -2,to = 2,length.out = K)
  net<-simulate_network(N = 10000,dyn = 0,dir = 0,theta = theta_true,pi = rep((1/K),K))
  
  #########################################################################################################
  ## Testing the EM undirected network code on the simulated dataset
  #undebug(wrapper_EM_undir)
  net_list<-wrapper_sc_SEM_undir(sim.net = net[[1]],nclust=2)
  
  cluster_ids_est<-as.vector(apply(X = matrix(1:1000),MARGIN = 1,FUN = function(x){
    cluster_id<-which.max(net_list[[1]][x,])
    return(cluster_id)
  }))
  length(which((net[[2]]==cluster_ids_est)==TRUE))
  
  BIC_results<-matrix(NA_real_,3,3)
  for (k in 2:4){
    net_list<-wrapper_EM_undir(sim.net = net[[1]],nclust=k,bandwidth = 4)
    BIC_results[k-1,]<-net_list[[2]]
  }
  param<-net_list[[1]]
  df_theta_est<-data.frame(1:51,param[[1]][[3]])
  names(df_theta_est)<-c("time","K1","K2")
  df_theta_est<-data.frame(melt(df_theta_est,id.vars = "time"),id="Estimated")
  df_theta_final<-rbind(df_theta,df_theta_est)
  p<-ggplot(data = df_theta_final)
  p1<-p+geom_point(mapping = aes(x = time,y = value,linetype=variable,color=id),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=variable,color=id))+labs(x="time",y="theta")
  return(p1)
}

#########################################################################################################
## Creating simulation module for simulating the network for directed case K>=2 and testing the EM_dir code
sim_module_EM_dir<-function(){
  ## Defining theta array T*K*2 directed case
  T_data<-51 ## Defining number of data points
  K<-2 ## Defining the number of clusters
  z <- seq(0,1,by=0.02)
  theta_true<-array(NA_real_,dim=c(T_data,K,2))
  for(i in 1:length(z)){
    theta_true[i,1,1]<-(1/2)*cos(3*pi*z[i])
    theta_true[i,2,1]<--((1/2)*cos(3*pi*z[i]))
    #theta_true[i,3,1]<--((1/2)*sin(3*pi*z[i]))-1.5
    theta_true[i,1,2]<-(1/2)*exp(z[i]) - 2
    theta_true[i,2,2]<-((1/2)*exp(1-z[i]) - 2)
    #theta_true[i,3,2]<-((1/2)*exp(z[i]) - 2)
  }
  # for(i in 1:length(z)){
  #   #theta_true[i,1,1]<-(1/2)*cos(3*pi*z[i])-0.75
  #   #theta_true[i,2,1]<--((1/2)*cos(3*pi*z[i])+1)-0.25
  #   #theta_true[i,2,1]<--((1/2)*sin(3*pi*z[i]))-1
  #   #theta_true[i,1,2]<-(1/2)*exp(z[i]) - 2.75
  #   #theta_true[i,2,2]<-((z[i]-0.5)^2)+((1/2)*exp(1/2)-2)+0.5
  #   #theta_true[i,2,2]<-((1/2)*exp(1-z[i]) - 2)
  # }
  #save(theta_true,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/theta_true/dir_K3_s1_theta_true.RData")
  
  df_theta_oe<-data.frame(1:51,index="oe",theta_true[,,1])
  df_theta_re<-data.frame(1:51,index="re",theta_true[,,2])
  df_theta<-rbind(df_theta_oe,df_theta_re)
  #names(df_theta)<-c("time","type","K1","K2","K3")
  names(df_theta)<-c("time","type","K1","K2")
  df_theta<-data.frame(melt(df_theta,id.vars = c("time","type")),id="true")
  p<-ggplot(data = df_theta)
  p_true<-p+geom_point(mapping = aes(x = time,y = value,linetype=type,color=variable),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=type,color=variable))+labs(x="time",y="theta")
  
  save(theta_true,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Simulation/dir_K2_s1_theta_true.RData"))
}
  #Simulating the network
  for (sim_index in 1:10){
    net_list<-simulate_network(N = 100,dyn = 1,dir = 1,theta = theta_true,pi = c(0.5,0.5))
    save(net_list,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Simulation/data/dir_K2_s1_n100_",sim_index+100,".RData"))
  }
  
  #ptm<-proc.time()
  #net_clust<-simulate_network(N = 100,dyn = 1,dir = 1,theta = theta_true,pi = c((1/3),(1/3),(1/3)))
  
  net_clust<-simulate_network(N = 100,dyn = 1,dir = 1,theta = theta_true,pi = c(0.5,0.5))
  net<-net_clust[[1]]
  cluster_ids_true<-net_clust[[2]]
  #########################################################################################################
  # debug(wrapper_EM_dir_CV)
  # result<-wrapper_EM_dir_CV(sim.net = net,nclust = 2, bandwidth = 4,fold_list = get_fold_list(sim.net = net,folds = 3),folds = 3)
  
  # Testing the directed network code on a simulated dataset
  net_list<-wrapper_EM_dir(sim.net = net,nclust=3,bandwidth = 4,sim_indicator = 1,theta_true = theta_true,K_true=2,cluster_ids_true=cluster_ids_true)
  
  # EM_dir_K3_s1_n100_h30<-list()
  # for (k in 1:4){
  #   EM_dir_K3_s1_n100_h30[[k]]<-wrapper_EM_dir(sim.net = net,nclust=k,bandwidth = 3.0,sim_indicator = 1,theta_true = theta_true,K_true=3,cluster_ids_true=cluster_ids_true)
  # }
  # save(EM_dir_K3_s1_n100_h30,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/results/EM_dir_K3_s1_n100_h30_",66,".RData"))
  # 
  # cluster_ids_est<-as.vector(apply(X = matrix(1:500),MARGIN = 1,FUN = function(x){
  #   cluster_id<-which.max(net_list[[1]][[1]][x,])
  #   return(cluster_id)
  # }))
  # length(which((true_clust_ids==cluster_ids)==TRUE))


  # BIC_results<-matrix(NA_real_,4,3)
  # for (k in 1:4){
  #   net_list<-wrapper_EM_dir(sim.net = net,nclust=k,bandwidth = 4,sim_indicator = 1,theta_true = theta_true,K_true=2,cluster_ids_true=cluster_ids_true)
  #   BIC_results[k,]<-net_list[[2]]
  # }
  
  param<-net_list[[1]]
  
  df_theta_oe_est<-data.frame(1:51,index="oe",param[[3]][,,1])
  df_theta_re_est<-data.frame(1:51,index="re",param[[3]][,,2])
  df_theta_est<-rbind(df_theta_oe_est,df_theta_re_est)
  names(df_theta_est)<-c("time","type","K2","K1")
  #names(df_theta_est)<-c("time","type","K1","K2")
  df_theta_est<-data.frame(melt(df_theta_est,id.vars = c("time","type")),id="Estimated")
  df_theta_final<-rbind(df_theta,df_theta_est)
  p_K1<-ggplot(data = df_theta_final[which(df_theta_final[,"variable"]=="K1"),])
  p_K2<-ggplot(data = df_theta_final[which(df_theta_final[,"variable"]=="K2"),])
  #p_K3<-ggplot(data = df_theta_final[which(df_theta_final[,"variable"]=="K3"),])
  p1<-p_K1+geom_point(mapping = aes(x = time,y = value,linetype=type,color=id),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=type,color=id))+labs(x="time",y="theta",title="K1")
  p2<-p_K2+geom_point(mapping = aes(x = time,y = value,linetype=type,color=id),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=type,color=id))+labs(x="time",y="theta",title="K2")
  #p3<-p_K3+geom_point(mapping = aes(x = time,y = value,linetype=type,color=id),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=type,color=id))+labs(x="time",y="theta",title="K3")
  # p<-ggplot(data = df_theta_est)
  p1
  p2
  p3
  # p+geom_point(mapping = aes(x = time,y = value,linetype=type,color=variable),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=type,color=variable))+labs(x="time",y="theta",title="K1")
  return(list(p1,p2))
}

#########################################################################################################
## Creating simulation module for simulating the network for directed case K>=2 and testing the EM_dir code
sim_module_EM_dir_K1<-function(){
  ## Defining theta array T*K*2 directed case
  T_data<-51 ## Defining number of data points
  K<-1 ## Defining the number of clusters
  z <- seq(0,1,by=0.02)
  theta_true<-array(NA_real_,dim=c(T_data,2))
  for(i in 1:length(z)){
    theta_true[i,1]<-1-sin(2*pi*z[i])
    theta_true[i,2]<-cos(3*pi*z[i])-1
  }
  
  df_theta_true<-data.frame("Time"=1:51,"Outgoing"=theta_true[,1],"Reciprocity"=theta_true[,2])
  df_theta_true<-data.frame(melt(df_theta_true,id.vars = "Time"),id="True")
  p<-ggplot(data = df_theta_true)
  p+geom_point(mapping = aes(x = Time,y = value,linetype=variable,color=variable),size=0.2)+geom_line(mapping = aes(x = Time,y = value,linetype=variable,color=variable))+labs(x="time",y="theta")
  
  save(theta_true,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/theta_true/dir_K1_s1_theta_true.RData")
  
  for (sim_index in 1:100){
    net_list<-simulate_network(N = 100,dyn = 1,dir = 1,theta = theta_true,pi = 1)
    save(net_list,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Simulation_new/data/dir/K1/s1/n100/dir_K1_s1_n100_",sim_index,".RData"))
  }
  
  ## Simulating the network
  net_clust<-simulate_network(N = 100,dyn = 1,dir = 1,theta = theta_true,pi = 1)
  net<-net_clust[[1]]
  cluster_ids_true<-net_clust[[2]]
  #########################################################################################################
  # Testing the directed network code on a simulated dataset
  net_list<-wrapper_EM_dir(sim.net = net,nclust=1,bandwidth = 2.5,sim_indicator = 1,theta_true = theta_true,K_true=1,cluster_ids_true=cluster_ids_true)
  
  param<-net_list[[1]]
  df_theta_est<-data.frame("Time"=1:51,"Outgoing"=param[,1],"Reciprocity"=param[,2])
  df_theta_est<-data.frame(melt(df_theta_est,id.vars = "Time"),id="Estimated")
  df_theta_final<-rbind(df_theta_est,df_theta_true)
  p<-ggplot(data = df_theta_final)
  p1<-p+geom_point(mapping = aes(x = Time,y = value,linetype=id,color=variable),size=0.2)+geom_line(mapping = aes(x = Time,y = value,linetype=id,color=variable))+labs(x="Time",y="Parameter",linetype="",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
  p1
  return(p1)
}

#########################################################################################################
## Creating simulation module for simulating the network for directed case K>=2 and testing the EM_dir code
sim_module_EM_dir_K1<-function(){
  ## Defining theta array T*K*2 directed case
  T_data<-51 ## Defining number of data points
  z <- seq(0,1,by=0.02)
  theta_true<-matrix(NA_real_,T_data,2)
  for(i in 1:length(z)){
    theta_true[i,1]<--(1/4)*cos(3*pi*z[i])-0.5
    #theta_true[i,2]<-(1/4)*cos(3*pi*z[i])-1.5
    theta_true[i,2]<-(1/2)*exp(z[i]) - 2
  }
  
  df_theta<-data.frame(1:51,theta_true)
  names(df_theta)<-c("time","oe","re")
  df_theta<-data.frame(melt(df_theta,id.vars = c("time")),id="true")
  p<-ggplot(data = df_theta)
  p_true<-p+geom_point(mapping = aes(x = time,y = value,color=variable),size=0.2)+geom_line(mapping = aes(x = time,y = value,color=variable))+labs(x="time",y="theta")
  
  ## Simulating the network
  net_clust<-simulate_network(N = 100,dyn = 1,dir = 1,theta = theta_true,pi = 1)
  net<-net_clust[[1]]
  cluster_ids_true<-net_clust[[2]]
  #########################################################################################################
  net_list<-wrapper_EM_dir(sim.net = net,nclust=2,bandwidth = 4,sim_indicator = 1,theta_true = theta_true,K_true=1,cluster_ids_true=cluster_ids_true)
  
  # cluster_ids_est<-as.vector(apply(X = matrix(1:100),MARGIN = 1,FUN = function(x){
  #   cluster_id<-which.max(net_list[[1]][[1]][x,])
  #   return(cluster_id)
  # }))
  # length(which((cluster_ids_true==cluster_ids_est)==TRUE))
  
  param<-net_list[[1]]
  df_theta_est<-data.frame(1:51,param)
  names(df_theta_est)<-c("time","oe","re")
  df_theta_est<-data.frame(melt(df_theta_est,id.vars = c("time")),id="Estimated")
  df_theta_final<-rbind(df_theta,df_theta_est)
  p<-ggplot(data = df_theta_final)
  p+geom_point(mapping = aes(x = time,y = value,linetype=variable,color=id),size=0.2)+geom_line(mapping = aes(x = time,y = value,linetype=variable,color=id))+labs(x="time",y="theta")
  return(list(p1,p2))
}

#########################################################################################################
#########################################################################################################
## International trade dataset
## Loading the internation trade dataset
load(file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Real Data/International trade data/Trade_data.Rdata")

####################################################
## Running the code over International trade dataset; parallelizing over both bandwidths and clusters
bandwidths<-c(2.5,3.0,3.5,4.0,4.5,5.0)
nclusts<-c(1,2,3,4,5,6)
band_clust_combn<-expand.grid(nclusts,bandwidths)
ITD_network_results_temp<-list()
cl <- makeCluster(6)
registerDoParallel(cl)
ITD_network_results_temp<-foreach (index = 1:nrow(band_clust_combn))%dopar% {
  wrapper_EM_undir(sim.net = trade50,nclust = band_clust_combn[index,1],bandwidth = band_clust_combn[index,2],sim_indicator = 0)
}
stopCluster(cl)

####################################################
## Compiling and saving the results
ITD_network_results<-list()
k<-1
for(i in 1:length(bandwidths)){
  ITD_network_results[[i]]<-list()
  for(j in 1:length(nclusts)){
    ITD_network_results[[i]][[j]]<-ITD_network_results_temp[[k]]
    k<-k+1
  }
}

save(ITD_network_results,file="/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Real Data/International trade data/ITD_network_results.RData")

####################################################
## NCV for International trade dataset
## Running the NCV code over International trade dataset; parallelizing over repeats, each repeat running several bandwidths
repeats<-6
bandwidths<-c(1.5,2,2.5,3,3.5,4)
band_chosen<-list()
cl <- makeCluster(6)
registerDoParallel(cl)
band_chosen<-foreach(repeat_index = 1:repeats)%dopar%{
  band_wrapper_EM_undir(sim.net = trade50,bandwidths = bandwidths,K=4,folds = 3)
}
stopCluster(cl)
bandwidths_chosen<-unlist(band_chosen)

####################################################
## Extracting the theta and gamma for chosen cluster and bandwidth
K_chosen<-4
bandwidth_chosen<-2.5
bandwidths<-c(2.5,3.0,3.5,4.0,4.5,5.0)
nclusts<-c(1,2,3,4,5,6)
load(file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Real Data/International trade data/ITD_network_results.RData")
result_final<-ITD_network_results[[which(bandwidths==bandwidth_chosen)]][[which(nclusts==K_chosen)]]
theta_final<-result_final[[1]][[3]]
gamma_final<-result_final[[1]][[1]]
pi_final<-result_final[[1]][[2]]
save(gamma_final,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Real Data/International trade data/gamma_final.RData")
save(pi_final,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Real Data/International trade data/pi_final.RData")

####################################################
## Plotting the theta
df_theta_final<-data.frame(1:20,theta_final)
names(df_theta_final)<-c("time","K1","K2","K3","K4")
p<-ggplot(data = df_theta_final)
p1<-p+geom_point(mapping = aes(x = time,y = K1))+geom_line(mapping = aes(x = time,y = K1))+labs(x="Time")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+coord_cartesian(ylim = c(-2.25, -1.85)) 
p2<-p+geom_point(mapping = aes(x = time,y = K2))+geom_line(mapping = aes(x = time,y = K2))+labs(x="Time")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+coord_cartesian(ylim = c(0.225, 0.625)) 
p3<-p+geom_point(mapping = aes(x = time,y = K3))+geom_line(mapping = aes(x = time,y = K3))+labs(x="Time")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+coord_cartesian(ylim = c(-0.78, -0.38)) 
p4<-p+geom_point(mapping = aes(x = time,y = K4))+geom_line(mapping = aes(x = time,y = K4))+labs(x="Time")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+coord_cartesian(ylim = c(2.05, 2.45)) 

multiplot(p1,p3,p2,p4,cols=2)

####################################################
## Plotting the clBIC checkmark plot
clBIC<-rep(NA_real_,length(nclusts))
for(k in 1:length(nclusts)){
  clBIC[k]<-ITD_network_results[[which(bandwidths==bandwidth_chosen)]][[k]][[2]][3]
}
df_clBIC<-data.frame("nclusts"=nclusts,"clBIC"=clBIC)
p<-ggplot(data = df_clBIC)
p+geom_point(mapping = aes(x = nclusts,y = clBIC))+geom_line(mapping = aes(x = nclusts,y = clBIC))+labs(x="Number of clusters",y="clBIC")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=15),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))

#########################################################################################################
#########################################################################################################
## Arms trade dataset
## Loading the arm trade dataset
load(file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/ArmTrade_raw _data.Rdata")

####################################################
## Preprocessing step
network_list<-list()
for (i in 1:20){
  network_list[[i]]<-eval(parse(text=paste0("net",i)))
}
node_list<-lapply(X = network_list,FUN = function(x){
  y<-colnames(x)
  return(y)
})
node_intersection<-Reduce(intersect,node_list)
node_intersection_len<-length(node_intersection)
net<-array(NA_real_,dim=c(node_intersection_len,node_intersection_len,20))
network_modified_list<-lapply(X = network_list,FUN = function(x){
  ids<-which(colnames(x)%in%node_intersection)
  y<-x[ids,ids]
  return(y)
})
net<-do.call(abind, list(network_modified_list,along=3))
str(net)

## Removing the unnecessary files
for (i in 1:20){
  eval(parse(text=paste0("remove(net",i,")")))
}

####################################################
## Screening out the isolated nodes

## Defining a function to remove the nodes that are isolated even in one time point from all time points
network_sieve<-function(df_list){
  ## Taking union of nodes for each dataframe
  node_union_list<-lapply(X = df_list,FUN = function(x){
    union_vector<-union(x$Node1,x$Node2)
    return(union_vector)
  })
  
  ## intersection of all vectors formed by individual unions. Every node in this vector must be present either in column 1 or column 2 of individual dataframes.
  node_labels<-node_intersection<-Reduce(intersect, node_union_list)
  
  ## subsetting each individual dataframe based on common nodes
  df_list_final<-lapply(X = df_list,FUN = function(x){
    ids_1<-c()
    ids_2<-c()
    for (i in 1:length(node_labels)){
      ids_1<-c(ids_1,which(x$Node1==node_labels[i]))
      ids_2<-c(ids_2,which(x$Node2==node_labels[i]))
    }
    ids<-intersect(ids_1,ids_2)
    y<-x[ids,]
    return(y)
  })
  
  # Creating network objects
  network_obj_list<-lapply(X = df_list_final, FUN = function(x){
    y<-data.frame(as.character(x[,1]),as.character(x[,2]))
    z<-as.network(y, directed=TRUE, matrix.type="edgelist")
    return(z)
  })
  
  ## Creating adjacency matrices
  network_adjmat_list<-lapply(X = network_obj_list,FUN = function(x){
    y<-as.matrix.network(x, matrix.type="adjacency")
    return(y)
  })
  
  network_adjmat_dim_list<-lapply(X = network_adjmat_list,FUN = function(x){
    y<-dim(x)
    return(y)
  })
  
  network_adjmat_dim_mat<-do.call(rbind,network_adjmat_dim_list)
  return(list(df_list_final,network_obj_list,network_adjmat_list,network_adjmat_dim_mat))
}

## Defining a sieve function to iterate multiple times the sieve function define above so that we get a network with common un-isolated nodes in all the time points.
network_sieve_iterator<-function(df){
  net_df<-network_sieve(df)
  index<-1
  while((length(unique(net_df[[4]][,1]))!=1)|(length(unique(net_df[[4]][,2]))!=1)){
    net_df<-network_sieve(net_df[[1]])
    index<-index+1
    print(index)
  }
  print(net_df[[4]])
  plot(net_df[[2]][[length(df)]])
  return(net_df)
}

df_net_edgelist<-list()
node_names<-attr(net,"dimnames")[[1]]
for (i in 1:20){
  df_net_edgelist[[i]]<-data.frame(as.matrix.network.edgelist(as.network(net[,,i], directed=TRUE, matrix.type="adjacency"),attrname = NULL))
  #df_net_edgelist[[i]]<-data.frame(as.edgelist(as.network(net[,,i], directed=TRUE, matrix.type="adjacency")))
  names(df_net_edgelist[[i]])<-c("Node1","Node2")
}
str(df_net_edgelist)

network_adjmat_list<-network_sieve_iterator(df_net_edgelist)[[3]]
node_size<-dim(network_adjmat_list[[1]])[1]
net<-array(NA_integer_,dim = c(node_size,node_size,length(network_adjmat_list)))
for (i in 1:length(network_adjmat_list)){
  net[,,i]<-network_adjmat_list[[i]]
}
str(net)

save(net,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/ArmTrade_un_isolated.Rdata")

####################################################
## Taking out all the isolated nodes reduces the node size a lot. Therefore we consider thresholding the number of un-isolated nodes to appear greater than or equal to threshold time points.
## Taking union of nodes for each dataframe
node_union_list<-lapply(X = df_net_edgelist,FUN = function(x){
  union_vector<-union(x$Node1,x$Node2)
  return(union_vector)
})
str(node_union_list)

node_union<-Reduce(union,node_union_list)
count<-rep(0,length(node_union))
for(i in 1:length(node_union)){
  for (j in 1:length(node_union_list)){
    if(node_union[i]%in%node_union_list[[j]]){
      count[i]<-count[i]+1
    }
  }
}

## Checking number of nodes for different thresholds
n_nodes<-rep(NA_integer_,20)
for (i in 1:20){
  n_nodes[i]<-length(which(count>=i))
}

list_updater<-function(df_list,node_labels){
  df_list_mod<-lapply(X = df_list,FUN = function(x){
    rownames(x)<-NULL
    total_nodes<-union(x$Node1,x$Node2)
    extra_nodes<-setdiff(total_nodes,node_labels)
    ## Removing the rows correspoding to extra nodes
    if (length(extra_nodes)>0){
      ids_1<-c()
      ids_2<-c()
      for (i in 1:length(extra_nodes)){
        ids_1<-c(ids_1,which(x$Node1==extra_nodes[i]))
        ids_2<-c(ids_2,which(x$Node2==extra_nodes[i]))
      }
      ids<-union(ids_1,ids_2)
      y<-x[-ids,]
      rownames(y)<-NULL
    }else{
      y<-x
    }
    return(y)
  })
  return(df_list_mod)
}
node_labels_updater<-function(df_list,node_labels,threshold){
  node_union_list<-lapply(X = df_list,FUN = function(x){
    union_vector<-union(x$Node1,x$Node2)
    return(union_vector)
  })
  
  node_union<-Reduce(union,node_union_list)
  count<-rep(0,length(node_union))
  for(i in 1:length(node_union)){
    for (j in 1:length(node_union_list)){
      if(node_union[i]%in%node_union_list[[j]]){
        count[i]<-count[i]+1
      }
    }
  }
  node_count_mat<-data.frame("NodeID"=node_union,"Count"=count)
  node_labels<-node_union[which(count>=threshold)]
  return(node_labels)
}
iterator<-function(df_list,node_labels,threshold){
  node_labels_prev<-0
  node_labels_curr<-node_labels
  index<-1
  while(identical(node_labels_prev,node_labels_curr)==FALSE){
    node_labels_prev<-node_labels_curr
    df_list<-list_updater(df_list = df_list,node_labels = node_labels)
    node_labels<-node_labels_updater(df_list = df_list,node_labels = node_labels,threshold = threshold)
    node_labels_curr<-sort(node_labels)
    index<-index+1
    #print(index)
  }
  ## Adding the missing nodes
  df_list_final<-lapply(X = df_list,FUN = function(x){
    total_nodes<-union(x$Node1,x$Node2)
    miss_nodes<-setdiff(node_labels,total_nodes)
    #print(length(miss_nodes))
    df_add<-data.frame(cbind(matrix(rep(miss_nodes,2),length(miss_nodes),2),rep(x$Year[1],length(miss_nodes))))
    names(df_add)<-c("Node1","Node2")
    z<-rbind(x,df_add)
    return(z)
  })
  return(list(df_list_final,node_labels))
}
## Loop over to thresholds to save different networks
for (thres_index in 3:20){
  threshold<-thres_index
  node_labels<-node_union[which(count>=threshold)]
  #debug(iterator)
  df_converged<-iterator(df_list = df_net_edgelist,node_labels = node_labels,threshold = threshold)
  df_list_final<-df_converged[[1]]
  # Creating network objects
  network_obj_list<-lapply(X = df_list_final, FUN = function(x){
    y<-data.frame(as.character(x$Node1),as.character(x$Node2))
    z<-as.network(y, directed=TRUE, matrix.type="edgelist")
    return(z)
  })

  ## Creating adjacency matrices
  network_adjmat_list<-lapply(X = network_obj_list,FUN = function(x){
    y<-as.matrix.network(x, matrix.type="adjacency")
    diag(y)<-0
    return(y)
  })
  
  node_size<-dim(network_adjmat_list[[2]])[1]
  net<-array(NA_integer_,dim = c(node_size,node_size,20))
  for (i in 1:20){
    net[,,i]<-network_adjmat_list[[i]]
  }
  node_names_curr<-node_names[as.numeric(rownames(network_adjmat_list[[1]]))]
  
  save(node_names_curr,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/threshold_data/ATD_node_names_",threshold,".RData"))
  
  save(net,file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/threshold_data/ATD_network_",threshold,".RData"))
  
  print(length(df_converged[[2]]))
}

#########################################################################################################
## Running the code over Arms trade dataset; parallelizing over both bandwidths and clusters
## Loading the threshold 10 dataset

for (thres_index in 3:20){
  load(file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/threshold_data/ATD_network_",thres_index,".RData"))
  bandwidths<-c(2.5,3.0,3.5,4.0,4.5,5.0)
  nclusts<-c(1,2,3,4,5,6)
  #temp<-wrapper_SEM_dir(sim.net = net,nclust = 5,bandwidth =3,sim_indicator = 0)
  band_clust_combn<-expand.grid(nclusts,bandwidths)
  ATD_network_results_temp<-list()
  cl <- makeCluster(6)
  registerDoParallel(cl)
  ATD_network_results_temp<-foreach (index = 1:nrow(band_clust_combn))%dopar% {
    wrapper_EM_dir(sim.net = net,nclust = band_clust_combn[index,1],bandwidth = band_clust_combn[index,2],sim_indicator = 0)
  }
  stopCluster(cl)
  
  ####################################################
  ## Compiling and saving the results
  ATD_network_results<-list()
  k<-1
  for(i in 1:length(bandwidths)){
    ATD_network_results[[i]]<-list()
    for(j in 1:length(nclusts)){
      ATD_network_results[[i]][[j]]<-ATD_network_results_temp[[k]]
      k<-k+1
    }
  }
  str(ATD_network_results)
  
  save(ATD_network_results,file=paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/results/ATD_network_results_un_isolated_threshold_",thres_index,".RData"))
}

bandwidths<-c(2.5,3.0,3.5,4.0,4.5,5.0)
nclusts<-c(1,2,3,4,5,6)
#temp<-wrapper_SEM_dir(sim.net = net,nclust = 5,bandwidth =3,sim_indicator = 0)
band_clust_combn<-expand.grid(nclusts,bandwidths)
ATD_network_results_temp<-list()
cl <- makeCluster(6)
registerDoParallel(cl)
ATD_network_results_temp<-foreach (index = 1:nrow(band_clust_combn))%dopar% {
  wrapper_EM_dir(sim.net = net,nclust = band_clust_combn[index,1],bandwidth = band_clust_combn[index,2],sim_indicator = 0)
}
stopCluster(cl)

####################################################
## Compiling and saving the results
ATD_network_results<-list()
k<-1
for(i in 1:length(bandwidths)){
  ATD_network_results[[i]]<-list()
  for(j in 1:length(nclusts)){
    ATD_network_results[[i]][[j]]<-ATD_network_results_temp[[k]]
    k<-k+1
  }
}
str(ATD_network_results)

save(ATD_network_results,file="/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/results/ATD_network_results_un_isolated_threshold_10.RData")

####################################################
## NCV for Arms trade dataset
## Running the NCV code over International trade dataset; parallelizing over repeats, each repeat running several bandwidths
load(file = paste0("/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/threshold_data/ATD_network_",3,".RData"))
str(net)
repeats<-6
bandwidths<-c(2,2.5,7.5,8,9,10)
band_chosen<-list()
cl <- makeCluster(6)
registerDoParallel(cl)
band_chosen<-foreach(repeat_index = 1:repeats)%dopar%{
  band_wrapper_EM_dir(sim.net = net,bandwidths = bandwidths,K=4,folds = 3)
}
stopCluster(cl)
bandwidths_chosen<-unlist(band_chosen)

####################################################
## Extracting the theta and gamma for chosen cluster and bandwidth
K_chosen<-4
bandwidth_chosen<-5
bandwidths<-c(2.5,3.0,3.5,4.0,4.5,5.0)
nclusts<-c(1,2,3,4,5,6)

#ATD_network_results<-wrapper_EM_dir(sim.net = net,nclust =4,bandwidth = 7.5,sim_indicator = 0)
load(file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/results/ATD_network_results_un_isolated_threshold_3.RData")

str(ATD_network_results)
result_final<-ATD_network_results[[which(bandwidths==bandwidth_chosen)]][[which(nclusts==K_chosen)]]
theta_final<-result_final[[1]][[3]]
gamma_final<-result_final[[1]][[1]]
clust_ids<-rep(NA_integer_,100)
for(i in 1:nrow(gamma_final)){
  clust_ids[i]<-which.max(gamma_final[i,])
}
sum(clust_ids==1)
sum(clust_ids==2)
sum(clust_ids==3)
sum(clust_ids==4)

round(gamma_final,2)
pi_final<-result_final[[1]][[2]]
save(gamma_final,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/results/threshold/3/bandwidth 5.0/gamma_final.RData")
save(pi_final,file = "/Users/Amal/Box Sync/PSU/Fall 2016/Research/Network Models/Project 1 (Semiparametric)/Real Data/Arms Trade Dataset/results/threshold/3/bandwidth 5.0/pi_final.RData")

####################################################
## Plotting the theta
#temp<-wrapper_EM_dir(sim.net = net,nclust = 3,bandwidth = 5,sim_indicator = 0)
df_theta_oe<-data.frame(1:20,ATD_network_results[[6]][[4]][[1]][[3]][,,1],id="oe")
df_theta_re<-data.frame(1:20,ATD_network_results[[6]][[4]][[1]][[3]][,,2],id="re")
#df_theta_oe<-data.frame(1:20,temp[[1]][[3]][,,1],id="Outgoing")
#df_theta_re<-data.frame(1:20,temp[[1]][[3]][,,2],id="Reciprocity")
df_theta_final<-rbind(df_theta_oe,df_theta_re)
names(df_theta_final)<-c("time","K1","K2","K3","K4","id")
p<-ggplot(data = df_theta_final)
p1<-p+geom_point(mapping = aes(x = time,y = K1,color=id))+geom_line(mapping = aes(x = time,y = K1,color=id))+labs(x="Time",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=10),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
p2<-p+geom_point(mapping = aes(x = time,y = K2,color=id))+geom_line(mapping = aes(x = time,y = K2,color=id))+labs(x="Time",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=10),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
p3<-p+geom_point(mapping = aes(x = time,y = K3,color=id))+geom_line(mapping = aes(x = time,y = K3,color=id))+labs(x="Time",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=10),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
p4<-p+geom_point(mapping = aes(x = time,y = K4,color=id))+geom_line(mapping = aes(x = time,y = K4,color=id))+labs(x="Time",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=10),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
#p5<-p+geom_point(mapping = aes(x = time,y = K5,color=id))+geom_line(mapping = aes(x = time,y = K5,color=id))+labs(x="Time",color="")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=10),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))

multiplot(p1,p3,p2,p4,cols=2)

####################################################
## Plotting the clBIC checkmark plot
bandwidth_chosen<-5
clBIC<-rep(NA_real_,length(nclusts))
for(k in 1:(length(nclusts))){
  clBIC[k]<-ATD_network_results[[which(bandwidths==bandwidth_chosen)]][[k]][[2]][3]
}
df_clBIC<-data.frame("nclusts"=nclusts[-1],"clBIC"=clBIC[-1])
p<-ggplot(data = df_clBIC)
p+geom_point(mapping = aes(x = nclusts,y = clBIC))+geom_line(mapping = aes(x = nclusts,y = clBIC))+labs(x="Number of clusters",y="clBIC")+theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),legend.text=element_text(size=10),strip.text.x = element_text(size = 20),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))



