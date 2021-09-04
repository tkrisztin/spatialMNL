# ### Bayesian logit model with Pólya Gamma prior ###
# ###  see Polson et al. (2013)

# 
# # install MASS package for multivariate normal random number
# # truncnorm for the truncated normal distribution
# # BayesLogit for Pólya Gamma random draws
 # if (!require("pacman")) install.packages("pacman")
 # pacman::p_load(MASS,truncnorm,BayesLogit)
 require(spdep)
 require(Matrix)
 require(MASS)
 require(truncnorm)
 require(BayesLogit)
 
   
######################################
#### DGP version
######################################

n=800
k = 7
p = 3

#X <- cbind( rnorm(n), rnorm(n),rnorm(n) )
#X <- cbind( runif(n),runif(n),runif(n) )
X <- matrix(rnorm(n*k),n,k)
X = scale(X,scale = F)

xy <- cbind(rnorm(n),rnorm(n))
nnb =  knn2nb(knearneigh(xy, k=6))
W=Matrix(nb2mat(nnb,style='W'))



# BETA = matrix(
#   c(sample(c(-5:5),prob = c(rep(.07,5),.3,rep(.07,5)),replace = T,size = k*(p-1)),
#     rep(0,k)),
#   k,p)

BETA = matrix(
  c(sample(c(-5:5),prob = c(rep(.07,5),.3,rep(.07,5)),replace = T,size = k*(p))),
  k,p)
BETA[,p] = 0

nn = rep(1,n) # ´number of trials

RHO = -0.5

AI = as.matrix(solve(diag(n) - RHO * W))
MU = AI %*% (X %*% BETA)
pr = exp(MU) / rowSums(exp(MU))
Y = pr
# Y = matrix(0,n,p)
# for (i in 1:n) {
#   pr2 =  sample(size=1,x = 3,prob = pr[i,])
#   Y[i,pr2] = 1
# }

vecAI = c(AI)
if (all(X[,1] == 1)) {
  ind_baseX = c(2:k)
} else {ind_baseX = c(1:k)}
DIRECT = matrix(0,length(ind_baseX),p)
INDIRECT = matrix(0,length(ind_baseX),p)
TOTAL = matrix(0,length(ind_baseX),p)
meanXs = apply(X,c(2),mean)
diag_ind = which(c(diag(n)) == 1)
ppp = 1
for (ppp in 1:p) {
 meanXmu = vecAI %*% t(BETA[ind_baseX,ppp] * meanXs[ind_baseX])
 meanXmu = exp(meanXmu)
 meanXmu = meanXmu / rowSums(meanXmu)

 bbb = vecAI %*% t(BETA[ind_baseX,ppp])
 pr_bbb = bbb
 for (kkk in 1:length(ind_baseX)) {
   pr_bbb[,kkk] = rowSums(meanXmu[,kkk] * vecAI %*% t(c(BETA[ind_baseX[kkk],])))
 }
 partial1 = meanXmu * (bbb - pr_bbb)

 DIRECT[,ppp] = colSums(partial1[diag_ind,])/n
 TOTAL[,ppp] = colSums(partial1) / n
 INDIRECT[,ppp] = TOTAL[,ppp] - DIRECT[,ppp]
}

######################################
#### DGP version  - done
######################################



n = nrow(X)
p = ncol(Y)
k = ncol(X)
 
nn = matrix(1,n,p)

#source("WmnlPGLogitMCMC.R")
#res = WmnlPGLogitMCMC(Y[,-p],matrix(1,n,1),X,1000,1000)

### let us assign prior values
# beta mean and variance are proper, but with high variance (= low prior information)
beta_prior_mean = matrix(0,k,p)
beta_prior_var = diag(k) * 10^4
beta_prob = function(rho,a) 1/beta(a,a) * ((1+rho)^(a-1) *(1-rho)^(a-1))/(2^(2*a - 1))
rho_a = 1.01


### set-up the gibbs sampler
# total number of draws
niter = 1000
# retain only S2 and discard S1, with S = S1 + S2
nretain = 500
ndiscard = niter - nretain
# save the posterior draws here
postb = array(0,c(k,p,nretain))
posty = array(0,c(n,p,nretain))
postom = array(0,c(n,p,nretain))
postr = matrix(0,1,nretain)

griddy_n = 100
source("lndetPaceBarry.R")
logdets = lndetPaceBarry(as.matrix(W),length.out = griddy_n+2)
logdets = logdets[-c(1,griddy_n + 2),]
rrhos = logdets[,2]
# Ais = array(0,c(n,n,griddy_n))
# AiXs = array(0,c(n,k,griddy_n))
# cat("Pre-calculate griddy GIBBS...")
# for (ii in 1:griddy_n) {
#   #Ais[,,ii] = as.matrix(solve(.sparseDiagonal(n) - rrhos[ii] * W))
#   tempA = .sparseDiagonal(n) - rrhos[ii] * W
#   Ai = solve(tempA)
#   AiXs[,,ii] = as.matrix(Ai %*% X)
#   #AiXKs[,ii] = t(AiXs[,,ii]) %*% kappa
#   #YAiXs[ii,] = t(Y) %*% AiXs[,,ii]
# }
# cat("Done!\n")
### calibration parameters for rho sampling
cc = 1 #scaling of rho proposals
c_adjust = 1.1 #proposal distribution adjustment
rho_accept = 0 #counter for rho acceptance rates
rtmp = matrix(0,1,niter) # rho acceptance rates


# starting values (won't matter after sufficient draws)
#curr.beta = mvrnorm(1,beta_prior_mean,beta_prior_var)
curr.beta = solve(crossprod(X)) %*% crossprod(X,Y)
curr.beta[,p] = 0
curr.xb = X %*% curr.beta
curr.om = matrix(0,n,p)
curr.rho = 0
curr.A = .sparseDiagonal(n) - curr.rho * W
curr.AiX = solve(.sparseDiagonal(n) - curr.rho * W) %*% X

# pre-calculate some terms for faster draws
beta_prior_var_inv = solve(beta_prior_var)
kappa = Y - nn/2

### Gibbs sampling
for (iter in 1:niter) {
  cat("iter:",iter,"curr.rho:",curr.rho,"\n")
  
  for (j in 1:(p-1)) {
    
    A = rowSums( exp( curr.AiX %*% curr.beta[,-j]) );
    c.j   = log(A);
    eta.j = as.matrix(curr.AiX %*% curr.beta[,j] - c.j);
    
    # sample omega
    curr.om[,j] = rpg(n, nn[,j], eta.j)
    
    # draw beta
    V = solve(beta_prior_var_inv + t(curr.AiX) %*%  (curr.AiX*curr.om[,j]) )
    b = V %*% (beta_prior_var_inv%*%beta_prior_mean[,j] +  t(curr.AiX) %*% (kappa[,j] + c.j * curr.om[,j]) )
    curr.beta[,j] = as.matrix(mvrnorm(1,b,V))
  }
  
  ## Metropolis-Hastings step for rho
  curr.llh = log(beta_prob(curr.rho,rho_a))
  for (j in 1:(p-1)) {
    A = rowSums( exp( curr.AiX %*% curr.beta[,-j]) );
    c.j   = log(A);
    ee1 = kappa[,j] / curr.om[,j] + c.j - curr.AiX %*% curr.beta[,j]
    curr.llh = curr.llh - sum(t(ee1 * curr.om[,j]) %*% ee1)/2
  }
  
  # curr.llh = log(beta_prob(curr.rho,rho_a))
  # for (j in 1:(p-1)) {
  #   A = rowSums( exp( curr.AiX %*% curr.beta[,-j]) );
  #   c.j   = log(A);
  #   eta.j = as.matrix(curr.AiX %*% curr.beta[,j] - c.j);
  #   exp_eta.j = exp(eta.j)
  #   ee1 = Y[,j] * (eta.j - log(1 + exp_eta.j)) + (nn[,j] - Y[,j]) * (1 - log(1 + exp_eta.j))
  #   curr.llh = curr.llh  + sum(ee1)
  # }

  accept = 0;
  while (accept!=1) {
    prop.rho = curr.rho + cc*rnorm(1,0,1)
    if (prop.rho<.95 && prop.rho>-.95) {
      accept = 1
    }
  }
  
  prop.A = .sparseDiagonal(n) - prop.rho * W
  prop.AiX = solve(prop.A) %*% X
  
  prop.llh = log(beta_prob(prop.rho,rho_a))
  for (j in 1:(p-1)) {
    A = rowSums( exp( prop.AiX %*% curr.beta[,-j]) );
    c.j   = log(A);
    ee1 = kappa[,j] / curr.om[,j] + c.j - prop.AiX %*% curr.beta[,j]
    prop.llh = prop.llh - sum(t(ee1 * curr.om[,j]) %*% ee1)/2
  }
  
  # prop.llh = log(beta_prob(prop.rho,rho_a))
  # for (j in 1:(p-1)) {
  #   A = rowSums( exp( prop.AiX %*% curr.beta[,-j]) );
  #   c.j   = log(A);
  #   eta.j = as.matrix(prop.AiX %*% curr.beta[,j] - c.j);
  #   exp_eta.j = exp(eta.j)
  #   ee1 = Y[,j] * (eta.j - log(1 + exp_eta.j)) + (nn[,j] - Y[,j]) * (1 - log(1 + exp_eta.j))
  #   prop.llh = prop.llh  + sum(ee1)
  # }
  
  acc_prob = min(1,exp(prop.llh - curr.llh))
  if (rbinom(1,1,acc_prob) == 1) {
    curr.rho = prop.rho
    rho_accept = rho_accept + 1
    #curr.logdet = prop.logdet
    curr.AiX = prop.AiX
  }
  # adjust candidate distribution based on acceptance probability
  rtmp[iter] = rho_accept/iter
  if (iter < ndiscard/2) {
    if (rtmp[iter]<.2) {
      cc <- cc/c_adjust
    } else if (rtmp[iter]>.4) {
      cc <- cc*c_adjust
    }
  }
  
  # we are past the burn-in, save the draws
  if (iter > ndiscard) {
    s = iter - ndiscard
    postb[,,s] = curr.beta
    curr.xb = X %*% curr.beta
    posty[,,s] = exp(curr.xb) / rowSums(exp(curr.xb))
    postom[,,s] = curr.om
    postr[s] = curr.rho
    
  }
}

### calculate posterior mean of beta and sigma
beta_mean_hat = apply(postb,c(1,2),mean)
y_mean_hat = apply(posty,c(1,2),mean)


thinning = 10
post_calcs = seq(1,nretain,by = thinning)
meanXs = apply(X,c(2),mean)
diag_ind = which(c(diag(n)) == 1)
if (all(X[,1] == 1)) {
  intercept = TRUE
} else {
  intercept = FALSE
}

if (intercept) {
  ind_baseX = c(2:k)
} else {
  ind_baseX = c(1:k)
}

ind_WX = ind_baseX + length(ind_baseX)

post.direct = array(0,c(length(ind_baseX),p,length(post_calcs)))
post.indirect = array(0,c(length(ind_baseX),p,length(post_calcs)))
post.total = array(0,c(length(ind_baseX),p,length(post_calcs)))

dimnames(post.direct)[[1]] = colnames(X)[ind_baseX]
dimnames(post.direct)[[2]] = colnames(Y)
dimnames(post.indirect)[[1]] = colnames(X)[ind_baseX]
dimnames(post.indirect)[[2]] = colnames(Y)
dimnames(post.total)[[1]] = colnames(X)[ind_baseX]
dimnames(post.total)[[2]] = colnames(Y)

# jjj = 1
# for (jjj in 1:length(post_calcs)) {
#   ind = post_calcs[jjj]
#   cat("Post direct/indirect - ",jjj,"of",length(post_calcs),"\n")
#   AI = as.matrix(solve(.sparseDiagonal(n) - postr[ind] * W))
#   vecAI = c(AI)
#   vecAIW = c(as.matrix(AI %*% W))
#   
#   ppp = 1
#   for (ppp in 1:p) {
#     meanXmu = vecAI %*% t(postb[ind_baseX,ppp,ind] * meanXs[ind_baseX]) #+ 
#       #vecAIW %*% t(postb[ind_WX,ppp,ind] * meanXs[ind_WX])
#     meanXmu = exp(meanXmu)
#     meanXmu = meanXmu / rowSums(meanXmu)
#     
#     bbb = vecAI %*% t(postb[ind_baseX,ppp,ind])  #+ vecAIW %*% t(postb[ind_WX,ppp,ind]) 
#     pr_bbb = bbb
#     for (kkk in 1:length(ind_baseX)) {
#       pr_bbb[,kkk] = rowSums(meanXmu[,kkk] * vecAI %*% t(c(postb[ind_baseX[kkk],,ind])))
#     } 
#     partial1 = meanXmu * (bbb - pr_bbb)
#   
#     post.direct[,ppp,jjj] = colSums(partial1[diag_ind,])/n
#     post.total[,ppp,jjj] = colSums(partial1) / n
#     post.indirect[,ppp,jjj] = post.total[,ppp,jjj] - post.direct[,ppp,jjj]
#   }
# }


