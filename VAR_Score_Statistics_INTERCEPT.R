## SCORE STATISTIC FOR VAR(p) MODEL


##EXAMPLES-------------------------
##example------------------------- bankfox data
bf_ts <- ts(bf[672:dim(bf)[1], c(1, 3:6, 8:10)])
bf_ts <- na.fill(bf_ts, fill = "extend")
ts.plot(bf_ts)

model <- ar(bf_ts, demean = T, method = "ols")
#Phi <- as.vector(t(as.matrix(model$ar)))
A_1 <-cbind(model$x.intercept,  matrix(model$ar, nrow=8, ncol=8))
eps <- model$resid

##example------------------------- change data
rData1 <- rSim(a1, e1)
rData1[1,] <- runif(5, 0, 0.02) #prevent NaN
rData2 <- rSim(a2, e2)
rData2[1,] <- runif(5, 0, 0.02)

var_change <- ts(rbind(rData1, rData2) + 3)
change_model <- ar(var_change, order.max = 1, demean = T, method = "ols")
a_change <- cbind(change_model$x.intercept, matrix( change_model$ar, nrow=5, ncol=5))
eps_change <- change_model$resid
plot(var_change)

##example------------------------- nochange data
nochange <- ts(rData1)
nochange_model <- ar(nochange, order.max = 1, demean = F, method = "ols")
a_nochange <- matrix(nochange_model$ar, nrow=5, ncol=5)
eps_nochange <- nochange_model$resid
plot(nochange)

##example------------------------- univariate data
univData1 <- arima.sim(n=2000, list(ar=c(a1_univ)), sd = 0.05) + 0.05; univData2 <- arima.sim(n=2000, list(ar=c(a2_univ)), sd = 0.05)
#univData1[1] <- univData2[1] <- 0.1 #prevent nan
univData <-  ts(c(univData1, univData2))
plot(univData)
univ_model <- ar(univData, aic=FALSE, order.max = 1, demean = T, method = "ols")
a_univ <- matrix(c(univ_model$x.intercept, univ_model$ar), nrow=1, ncol=2)
eps_univ <- univ_model$resid





########################################
##SCORE----------------------------------

##estimating fn for channel i, time k
getH_ik <- function(x, i,k, Phi, eps) { 
  if(is.null(dim(x))) {x <- matrix(x); Phi <- matrix(Phi); eps <- matrix(eps)} #handle univariate case
  a <- t(Phi[i,]) #A_1[i,]
  X <- c(1, x[(k-1),]) #with intercept
  #X <- t(X)
  y <- as.numeric(x[k,i])
  e <- eps[k,i]
  H_ik <- t(- y*X + a%*%X%*%t(X) - e*X)
  return( H_ik)
  #return(list(X%*%t(X), cov(X,X)))
}
#getH_ik(bf_ts, i=1, k=10, Phi = A_1, eps = (eps))
#getH_ik(matrix(univData), 1, 100, (a_univ), matrix(eps_univ))

##make whole estimating function H (pd^2 vector) at time k
makeH_k <- function(x, k, p, Phi, eps){ 
  d <- dim(x)[2]
  if(is.null(dim(x))) {d <- 1} #handle univariate case
  H <- t(rep(0, d + d*d*p)) #accounts for intercept
  for (i in 1:d) {
    H[((i-1)*(d+1)+1):((i-1)*(d+1)+ d+1)] <- getH_ik(x,i,k,Phi,eps) #accounts for intercept
  }
  return(t(H))
}
#makeH_k(bf_ts, k=10, p=1, Phi = A_1, eps = eps)
#makeH_k(matrix(univData), k=101, p=1, Phi = (a_univ), eps = matrix(eps_univ))


##evaluate estimating function H at all time steps
makeH_all <- function(x, p, G, Phi, eps){ 
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  if(is.null(dim(x))) {n <- length(x); d<- 1} #handle univariate case
  K <- (G+1):(n-G) #time indices to evaluate over
  H_all <- matrix(0, nrow = d+d*d*p, ncol = n) #matrix of H values #accounts for intercept
  for (t in K) {
    H_all[,t] <- makeH_k(x, k=t, p, Phi, eps) 
  }
  return(H_all)
}
#H_all_bfts <- makeH_all(x=bf_ts, p=1, G=100, Phi=A_1, eps=eps)
#H_all_nochange <- makeH_all(x=nochange, p=1, G=10, Phi=a_nochange, eps=eps_nochange)
H_all_univ <- makeH_all(matrix(univData), 1, 100, (a_univ), matrix(eps_univ))

get_DiagH_rootinv <- function(x, k, G, H_all){
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  l <- (k-G):(k-1) #lower time index
  u <- k:(k+G-1) #upper time index
  H_list <- list(1: d )
  for (i in 1:d) {
    #if(d>1){ #d>1
      H_u <- H_all[((i-1)*(d+1)+1):((i-1)*(d+1)+ d+1),u] #accounts for intercept
      H_l <- H_all[((i-1)*(d+1)+1):((i-1)*(d+1)+ d+1), l] 
      Hbar_u <- colMeans(H_u) 
      H_u_centred <- (H_u - Hbar_u)
      Hbar_l <- colMeans(H_l)
      H_l_centred <- (H_l - Hbar_l)
      H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred) #sum  } else Hbar_u <- matrix(mean(H_u))
    # } else { #d=1
    #   H_u <- H_all[,u]
    #   H_l <- H_all[, l] 
    #   Hbar_u <- mean(H_u)
    #   H_u_centred <- (H_u - Hbar_u) 
    #   Hbar_l <- mean(H_l)
    #   H_l_centred <- (H_l - Hbar_l)
    #   H_out <- t(H_l_centred) %*% (H_l_centred) + t(H_u_centred) %*% (H_u_centred)
    # }
    ##eigen decomposition
    e <- eigen(H_out) #SVD of H_out
    V <- e$vectors
    if(d>1) {
      H_ <- V %*% diag((e$values)^(-0.5)) %*% t(V)
    } else {H_ <- e$values^(-0.5) * V %*% t(V)}
    H_list[[i]] <- H_
  }
  Sig_ <-  sqrt(2*G) *Matrix::bdiag(H_list) #coerce into block diagonal form
  return(Sig_)
}
#get_DiagH_rootinv(x=as.matrix(bf_ts), k =200,G=100, H_all=H_all_bfts)
#get_DiagH_rootinv(x=as.matrix(univData), k =200,G=100, H_all= (H_all_univ))

get_FullH_rootinv <- function(x, k, G, H_all){
  if(is.null(dim(x))) {x <- matrix(x); H_all = matrix(H_all)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  l <- (k-G):(k-1) #lower time index
  u <- k:(k+G-1) #upper time index
#    if(d>1){ #d>1
      H_u <- H_all[,u]
      H_l <- H_all[,l] 
      Hbar_u <- colMeans(H_u) 
      H_u_centred <- (H_u - Hbar_u)
      Hbar_l <- colMeans(H_l)
      H_l_centred <- (H_l - Hbar_l)
      H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred) #sum  } else Hbar_u <- matrix(mean(H_u))
    # } else { #d=1
    #   H_u <- H_all[u]
    #   H_l <- H_all[l] 
    #   Hbar_u <- mean(H_u)
    #   H_u_centred <- (H_u - Hbar_u) 
    #   Hbar_l <- mean(H_l)
    #   H_l_centred <- (H_l - Hbar_l)
    #   H_out <- t(H_l_centred) %*% (H_l_centred) + t(H_u_centred) %*% (H_u_centred)
    # }
    ##eigen decomposition
    e <- eigen(H_out) #SVD of H_out
    V <- e$vectors
    if(d>1) {
      H_ <- V %*% diag((e$values)^(-0.5)) %*% t(V)
    } else {H_ <- e$values^(-0.5) * V %*% t(V)}
  Sig_ <-  sqrt(2*G) * H_ 
  return(Sig_)
}
#get_FullH_rootinv(x=as.matrix(bf_ts), k =200,G=100, H_all=H_all_bfts)
#get_FullH_rootinv(x=as.matrix(univData), k =200,G=100, H_all= (H_all_univ))

##Difference Vector at k
getA <- function(x, k, p, G, Phi, eps, H_all){ 
  d <- dim(x)[2]
  r <- H_all[,(k+1):(k+G)]; l <- H_all[,(k-G+1):(k)] #left/right of window
  #if (d > 1){ 
    A <- rowSums(r) - rowSums(l)#difference 
  #} else {A <- sum(r) - sum(l)}
  A <- as.matrix(A)
  return(A)
}
#getA(x=bf_ts, k = 200, p=1, G=100, Phi=A_1, eps=eps, H_all = H_all_bfts)
#getA(x=nochange, k = 10, p=1, G=10, Phi=a_nochange, eps=eps_nochange, H_all = H_all_example)
#getA(matrix(univData), 200, 1, 100, (a_univ), (eps_univ), H_all_univ)

##sigma_i options---------------------------------------------
##global estimate for sigma^2_i 
getsigma_iGlobal <- function(eps,i){
  if(is.null(dim(eps))) {eps <- matrix(eps)} #handle univariate case
  mean(eps[-1,i]^2)
}
#getsigma_iGlobal(eps,1)
#getsigma_iGlobal(eps_univ,1)

##global estimate, all i=1,...,d, for sigma^2
getsigma_dGlobal <- function(eps){
  d <- dim(eps)[2]
  if(is.null(dim(eps))) d<- 1 #handle univariate case
  sigma_d <- rep(0,d)
  for (i in 1:d) { #sigma(i)
  sigma_d[i] <- getsigma_iGlobal(eps,i)
  } 
  return(sigma_d)
}
#bf_sigma_d <- getsigma_dGlobal(eps)
#getsigma_dGlobal(eps_univ)


##inverse root of global variance as block-diag matrix, altenative to sigma_root_inv  ####NOT USED - SEE NEXT FN
# get_S_ <- function(x, eps, p){ 
#   n <- dim(x)[1]
#   d <- dim(x)[2]
#   sigma_d <- rep(0,d)
#   S_list <- list(1:d)
#   for (i in 1:d) { 
#     sigma_d[i] <- getsigma_iGlobal(eps,i)^(-0.5)#sigma(i)^{-1/2}
#     S_list[[i]] <- diag(sigma_d[i], d,d) #as diagonal matrix
#   }
#   S_ <-  Matrix::bdiag(S_list) #coerce into block diagonal form
#   return(S_)
# }
#S_example <- get_S_(x= bf_ts, eps = eps,p=1)

##estimate Sigma_k^{-1/2} inverse of long-run covariance at k with DiagC
get_DiagC_rootinv <- function(x, eps, sigma_d, k, G){ 
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  xk <-  cbind(1, x[(k-G):(k+G-1),]) #lower time sample #includes intercept
  C <- 1/(2*G) * t(xk) %*% (xk)
  ##eigen decomposition
  e <- eigen(C) #SVD of C
  V <- e$vectors
  #if(d>1) {
    C_ <- V %*% diag( e$values^(-0.5) ) %*% t(V)
  #} else {C_ <- e$values^(-0.5) * V %*% t(V)}
  C_list <- list(1:d)
  for (i in 1:d) {
    C_list[[i]] <- sigma_d[i]^(-0.5) * C_
  }
  Sig_ <-  (2)^(-0.5) *Matrix::bdiag(C_list) #coerce into block diagonal form
  return(Sig_)
}
#Sig_example <- get_DiagC_rootinv(as.matrix(bf_ts), eps, bf_sigma_d, 400,100)
#get_DiagC_rootinv(matrix(univData), matrix(eps_univ), getsigma_dGlobal(eps_univ), 400, 100)


## inverse sqrt of local uncentred covariance C_nk
# get_inv_root_C_nk <- function(x, p, l, u){
#   n <- dim(x)[1]
#   d <- dim(x)[2]
#   xk <- x[(l-1):(u-1),] #time sample
#   c <- 1/(u-l) * t(xk) %*% (xk) #uncentred covariance of sample
#   e <- eigen(c) #SVD of c
#   V <- e$vectors
#   c_ <- V %*% diag((e$values)^(-.5)) %*% t(V) #c^{1/2}
#   c_list <- list(1:d)
#   for (i in 1:d) {
#     c_list[[i]] <- c_
#   }
#   C_ <-  Matrix::bdiag(c_list) #coerce into block diagonal form
#   return(C_)
# }
#get_inv_root_C_nk(bf_ts, p=1, l= 10, u= 100)

##evaluate statistic at time k
getTkn <- function(x, k, p, G, Phi, eps, H_all, sigma_d, estim = "DiagH"){ 
  A <- getA(x, k, p, G, Phi, eps, H_all) #difference matrix
  ##Sigma estimator options------
  if(estim == "DiagC") Sig_ <- get_DiagC_rootinv(x,eps,sigma_d,k,G)
  if(estim == "DiagH")Sig_ <- get_DiagH_rootinv(x,k,G,H_all) #DiagH estimator for Sigma
  if(estim == "FullH")Sig_ <- get_FullH_rootinv(x,k,G,H_all) #FullH estimator
  #------------------------------
  Tkn <- (2*G)^(-.5) * as.numeric( sqrt(t(A) %*% Sig_ %*% A) )#* norm(T_in,type = "F")
  return(Tkn)
}
#getTkn(x=as.matrix(bf_ts), k=200, p=1, G=100, Phi=A_1, eps=eps, H_all = H_all_bfts, sigma_d = bf_sigma_d) 
#getTkn(matrix(univData), k=200, p=1, G=100, Phi=a_univ, eps=eps_univ, H_all=)

##get statistic over all times K
getT <- function(x, p, G, Phi, eps, estim){ 
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  K <- (G+1):(n-G) #time indices to evaluate over
  H_all <- makeH_all(x, p, G, Phi, eps)
  Tkn <- rep(0,n)
  sigma_d <- getsigma_dGlobal(eps)
   #Inverse root long-run covariance
  #S_ <- get_S_(x,eps,p)
  for(k in K){
    Tkn[k] <- getTkn(x,k,p,G,Phi,eps, H_all, sigma_d, estim = estim)#Sig_) 
  }
  return(Tkn)
  #return(Sig_)
}
#Tn_example <- getT(x=bf_ts, p=1, G= 100, Phi=A_1, eps=eps, estim="DiagH")
#getT(x=var_change, p=1, G=40, Phi = a_change, eps = eps_change)

##get change point estimates
get_cps <- function(Tn, D_n, G, nu = 1/4){
    n <- length(Tn)
    rshift <- c(Tn[-1],0); lshift <- c(0,Tn[-n]); 
    over <- (Tn >D_n) #indices are greater than D_n?
    v <- which(over & lshift < D_n) #lowers
    #v <- which(v) #append 0
    w <- which(over & rshift < D_n) #uppers
    #w <- which(w) #append n
    nu_remove <- which(w-v >= nu*G) #nu(epsilon) test for distance between 
    v_nu <- v[nu_remove]; w_nu <- c(w[nu_remove],n)
    q <- length(v_nu) #number of CPs
    cps <- rep(0, q)
    for (i in 1:q) {
      cps[i] <- v_nu[i] + which.max(Tn[ (v_nu[i]):(w_nu[i]) ] )
    }
    return(cps)
}
get_cps(Tn = Tn_example, D_n = 5, G= 100,nu = 1/4)

##Score-type test
test_Score <- function(x, p, G, Phi, eps, alpha = 0.05, estim="DiagH"){ 
  if(is.null(dim(x)) || dim(x)[2] == 1 ) {x <- matrix(x); Phi <- matrix(Phi); eps <- matrix(eps)} #handle univariate case
  n <- dim(x)[1] #dimensions
  d <- dim(x)[2] 
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d/2 * log(log(n/G)) - log(2/3 * gamma(d/2))
  D_n <- (b+c_alpha)/a #threshold
  Reject <- FALSE
  ##Run test-----------------------------
  Tn <- ts(getT(x,p,G,Phi,eps,estim)) #evaluate statistic at each time k
  test_stat <- max(Tn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Tn,D_n,G, nu=1/4)
  } 
  ##Plot------------------------------------
  plot(Tn) # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = Tn, cps = cps, plot = pl, estim = estim)
  return(out)
}


##  bf example
bf_test <- test_Score(x=bf_ts, p=1, G=200, Phi = A_1, eps = eps, alpha = 0.05)
bf_test
##  change example
change_test <- test_Score(x=var_change, p=1, G=200, Phi = a_change, eps = eps_change, alpha = 0.05, estim = "DiagC") 
change_test$plot
##  nochange example
nochange_test <- test_Score(x=nochange, p=1, G=150, Phi = a_nochange, eps = eps_nochange, alpha = 0.1) 
nochange_test
## univariate example
univ_test <- test_Score(x=matrix(univData), p=1, G= 300, Phi = (a_univ), eps= matrix(eps_univ), alpha=0.1)
univ_test

