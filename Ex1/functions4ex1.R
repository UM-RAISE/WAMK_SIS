library(mvtnorm)
library(cubature)
# function to get X
get_X = function(num,mu_lst,sigma_lst,pi){
  K = length(pi)
  components = sample(1:K, prob = pi, size = num, replace = T)
  sample_X = c()
  for (n in 1:K){
    num_K = sum(components==n)
    if (num_K==0){
      next
    }
    sample_X = rbind(sample_X,rmvnorm(num_K, mean = mu_lst[[n]], sigma = sigma_lst[[n]]))
  }
  return(sample_X)
}
# fucntion to get Y
mu = function(X){
  X1 = X[,1]
  X2 = X[,2]
  X3 = X[,3]
  m = 40*(1-exp(-0.2*sqrt((X1^2+X2^2)/2)))+
    20*(1-exp(-0.2*abs(X1)))+
    5*(1-exp(-0.2*sqrt((X2^2+X3^2)/2)))-
    exp(cos(2*pi*X1*X3))-
    exp(cos(2*pi*X2*X3))-
    exp(cos(2*pi*X1*X2))-
    exp(cos(2*pi*X1*X2*X3))
  # exp((cos(2*pi*X1)+cos(2*pi*X3))/2)-
  # exp((cos(2*pi*X2)+cos(2*pi*X3))/2)
  return(m)
}

get_Y = function(X){
  n = nrow(X)
  Y = rnorm(n,mu(X),1)
  return(Y)
}



# sampling X using acceptance-rejection algorithm
sampling_X = function(num,X,Z,h_mat,combo,weight){
  D = ncol(X)
  x = c()
  nrow(x)
  i=0
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    x_new = rmvnorm(1,rep(0,D))
    u = runif(1,min=0,max=dmvnorm(x_new))
    if (u<=sqrt(S_add(x_new,X,Z,h_mat,combo,weight))*dmvnorm(x_new)){
      x = rbind(x,x_new)
    }
  }
  return(x)
}



get_Est = function(Y,H,l){
  prob = sum((Y>l)*1/H/sum(1/H))
  return(prob)
}
# # calculate C
# integrand_C = function(x,X_mat,Z,h1,h2){
#   x_mat = matrix(x,ncol = 2)
#   int = sqrt(S_bi(x_mat,X_mat,Z,h1,h2))*dmvnorm(x_mat)
#   return(int)
# }
# 
# get_C = function(X_mat,Z,h1,h2){
#   C = (adaptIntegrate(function(x) integrand_C(x,X_mat,Z,h1,h2),
#                       lowerLimit = c(-Inf,-Inf),
#                       upperLimit = c(Inf,Inf)))$integral
#   return(C)
# }





