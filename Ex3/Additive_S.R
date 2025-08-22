
source('AMISE_bivariate.R')

S_add = function(x_mat, X_mat, Z, h_mat, combo,weight){
  S_est = 0
  for (c in 1:nrow(combo)){
    cb = combo[c,]
    h = h_mat[c,]
    w = weight[c]
    S_temp = S_bi(matrix(x_mat[,cb],ncol=2),X_mat[,cb],Z,h[1],h[2])
    S_est = S_est+w*S_temp
  }
  S = S_est
  S[is.na(S)] = 0
  S[S==Inf] = 0
  return(S)
}