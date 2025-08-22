library(doParallel)
library(gplots)
######### This is Kernel SIS2
rm(list = ls())
set.seed(777)
source('AMISE_bivariate.R')
source('functions4ex1.R')
source('Additive_S.R')
#------ set up Parallel
ncores=2
cl = makeCluster(ncores)

registerDoParallel(cl)



#------ start experiments
POE_est = foreach (j = 1:100, .combine = c, .errorhandling = 'remove') %dopar%{
  set.seed(j+123)
  library(pracma)
  library(mvtnorm)
  library(cubature)
  # parameters initialization
  D = 3 # Dimension
  combo = t(combn(seq(1:D),2))
  Nm = 50
  Nm_init = 50
  l_target = 17.8979 #8.95 for D=3, 8.70 for D=4
  l = l_target
  iter = 5
  # mu_0 = list(rep(0,D))
  # sigma_0 = list(diag(4,D))
  # pi_0 = 1
  # Get initial samples
  sample_X = matrix(runif(Nm*D,-5,5),ncol = D)
  sample_Y = get_Y(sample_X)
  sample_Z = (sample_Y>l)
  sample_X_all = sample_X
  sample_Z_all = sample_Z
  # get starting points of bandwidth
  h_init_vec = c()
  for (i in 1:D){
    h = (optimize(function(x) AMISE(sample_X[,i],sample_Z,x,-5,5),interval = c(0.1,1.5)))$minimum
    h_init_vec[i] = h
  }
  h_init_mat = matrix(h_init_vec[combo],ncol = 2)
  
  est_vec = c()
  for (simu in 1:5){
    # Choose optimal bandwidth
    h_new_mat = 0.8*h_init_mat
    for (c in 1:nrow(combo)){
      cb = combo[c,] 
      h_init = h_init_mat[c,]
      h = try(get_h_V2_final(sample_X_all[,cb],sample_Z_all,h_init[1],h_init[2]),silent = T)
      if("try-error" %in% class(h)){
        next
      }else{
        h_new = h[nrow(h),]
        h_new_mat[c,] = h_new
      }
    }
    # deciding weights
    cost_vec = c()
    for (c in 1:nrow(combo)){
      cb = combo[c,] 
      h = h_new_mat[c,]
      S_temp = S_bi(sample_X_all[,cb],sample_X_all[,cb],sample_Z_all,h[1],h[2])
      cost = -sum(sample_Z_all*log(S_temp)+(1-sample_Z_all)*log(1-S_temp), na.rm = T)
      cost_vec[c] = cost
    }
    weight = 1/cost_vec/sum(1/cost_vec)
    # get new samples for the next round
    sample_X = sampling_X(Nm,sample_X_all,sample_Z_all,h_new_mat,combo,weight)
    #hist2d(sample_X, same.scale = T)
    sample_Y = get_Y(sample_X)
    sample_Z = (sample_Y>l)
    sample_S = S_add(sample_X,sample_X_all,sample_Z_all,h_new_mat,combo,weight)
    sample_H = sqrt(sample_S)
    est_temp = get_Est(sample_Y,sample_H,l_target)
    
    est_vec = c(est_vec,est_temp)
    # update sample bank
    sample_X_all = rbind(sample_X_all,sample_X)
    sample_Z_all = c(sample_Z_all,sample_Z)
    h_init_mat = h_new_mat
  }
  prob = mean(est_vec)
  prob
}

stopCluster(cl)
