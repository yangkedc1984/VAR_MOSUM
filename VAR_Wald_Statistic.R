## WALD STATISTIC FOR VAR(p) MODEL



##EXAMPLES-------------------------



########################################
##WALD-------------------------------

##local (l,u) regression parameter for channel i
get_a_lu_i <- function(x, i, l, u){
  y <- x[l:u,i] #response
  Xt <- t(x[(l-1):(u-1),]) #regressors ##p = 1 case only
  y_soln <- colSums(diag(y) %*% t(Xt))
  X_soln <- solve(Xt %*% t(Xt))
  a_lu_i <- t(y_soln %*% X_soln) 
  return(a_lu_i)
} 
get_a_lu_i(bf_ts, 1, 10, 100)

##local (l,u) regression parameter for all channels
make_a_lu <- function(x, p, l, u){
  d <- dim(x)[2]
  a_lu <- rep(0, d*d*p)
  for (i in 1:d) {
    a_lu[((i-1)*d+1):( (i-1)*d+d)] <- get_a_lu_i(x,i,l,u)
  }
  return(a_lu)
}
make_a_lu(bf_ts, p=1, l= 10, u= 100)

## sqrt of local uncentred covariance C_nk
get_V_nk <- function(x, p, l, u){
  d <- dim(x)[2]
  xk <- x[(l-1):(u-1),] #time sample
  C <- 1/(u-l) * t(xk) %*% (xk) #uncentred covariance of sample
  Vlist <- list(1:d)
  for(i in 1:d) {
    Vlist[[i]] <- 2 *(C) #coerce into block diagonal form
  }
  out <- Matrix::bdiag(Vlist)
  return(out)
}
#get_V_nk(as.matrix(bf_ts), p=1, l= 10, u= 100)

## LOCAL1 estimator of variance in channel i at k
getsigma_i_kLOCAL1 <- function(x, i, k, G, a_upper, a_lower) {
  x_upper <- t(x[(k):(k+G-1),]) #upper sample
  res_upper <-  x[(k+1):(k+G), i] - t(a_upper) %*% x_upper #upper residuals
  x_lower <- t(x[(k-G):(k-1),]) #lower sample
  res_lower <-  x[(k-G+1):(k), i] - t(a_lower) %*% x_lower #lower residuals
  sigma_i <- 1/(2*G) * (mean(res_upper^2) +mean(res_lower^2) ) #LOCAL1
  return(sigma_i)
}
#a_upper_example <- get_a_lu_i(bf_ts, i=1, l= 100, u= 130)
#a_lower_example <- get_a_lu_i(bf_ts, i=1, l= 69, u= 99)
#getsigma_i_kLOCAL1(x = bf_ts, i=1, k = 100, G= 30, a_upper = a_upper_example, a_lower = a_lower_example)


## sigma^2_i estimate for all times k; example function only, not used in test
getsigma_i_all <- function(x,i,G){
  n <- dim(x)[1]
  d <- dim(x)[2]
  K <- (G+3):(n-G-2)
  out <- rep(0, n)
  for (k in K) {
    a_i_upper <- get_a_lu_i(x, i, l=k+1, u=k+G)
    a_i_lower <- get_a_lu_i(x, i, l=k-G+1, u=k)
    out[k] <- getsigma_i_kLOCAL1(x,i,k,G, a_i_upper, a_i_lower)
  }
  return(ts(out))
}
#plot(1/sqrt(getsigma_i_all(x=bf_ts,i=1,G=50)) )
#plot(1/sqrt(getsigma_i_all(x=var_change,i=1,G=50)) ); abline( h= getsigma_iGlobal(eps,1)^(-0.5) )


#get_Gamma <- function(){
#}

## evaluate Wkn (test statistic) for time k
get_Wkn <- function(x, p, k, G, H_all, estim){ 
  d <- dim(x)[2]
  V <- get_V_nk(x, p, k-G+1, k)
    
  a_upper <- make_a_lu(x, p, l=k+1, u=k+G)
  a_lower <- make_a_lu(x, p, l=k-G+1, u=k)
  ##Sigma estimator options------
  #if(estim == "DiagC")Sig_ <- get_DiagC_rootinv(x,eps,sigma_d,k,G)
  if(estim == "DiagH")Sig_ <- get_DiagH_rootinv(x,k,G,H_all) #DiagH estimator for Sigma
  if(estim == "FullH")Sig_ <- get_FullH_rootinv(x,k,G,H_all) #FullH estimator
  #------------------------------
  #W <- W + 1/sigma_i * norm(root_C %*% (a_i_upper - a_i_lower), type = "F")
  W_mat <- as.matrix(Sig_ %*%V %*% (a_upper-a_lower)) #argument for norm
  W <- sqrt(G/2) * norm(W_mat, type="F")
  return(W)
}
get_Wkn(x = as.matrix(bf_ts), p=1, k= 200, G=100, H_all=H_all_bfts, estim="DiagH")

get_W <- function(x, p, G, Phi, eps, estim){
  n <- dim(x)[1]
  K <- (G+1):(n-G)
  H_all <- makeH_all(x, p, G, Phi, eps)
  out <- rep(0, n)
  for (k in K) {
    out[k] <- get_Wkn(x,p,k,G,H_all,estim)
  }
  return(out)
}
get_W(x= as.matrix(bf_ts), p=1, G=100, Phi=A_1, eps=eps, estim="DiagH")

##get change point estimates -- same as Score procedure


##Wald-type test
test_Wald <- function(x, p, G, Phi, eps, alpha = 0.05, estim="DiagH"){ 
  if(is.null(dim(x)) || dim(x)[2] == 1 ) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1] #dimensions
  d <- dim(x)[2] 
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d/2 * log(log(n/G)) - log(2/3 * gamma(d/2))
  D_n <- (b+c_alpha)/a #threshold
  Reject <- FALSE
  ##Run test-----------------------------
  Wn <- ts(get_W(x,p,G,Phi,eps,estim)) #evaluate statistic at each time k
  test_stat <- max(Wn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Wn,D_n,G, nu=1/4)
  } 
  ##Plot------------------------------------
  plot(Wn) # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = Wn, cps = cps, plot = pl, estim=estim)
  return(out)
}

## examples
# bf
bf_test_wald <- test_Wald(x= as.matrix(bf_ts), p=1, G=200, Phi = A_1, eps=eps, alpha= .1, estim = "DiagH")
bf_test_wald
# change
change_test_wald <- test_Wald(x= var_change, p=1, G= 200, Phi = a_change, eps=eps_change, alpha = .1, estim = "DiagH")
change_test_wald
# no change
nochange_test <- test_Wald(x=nochange, p=1, G=150, Phi = a_nochange, eps = eps_nochange, alpha = 0.1) 
nochange_test
#univariate
univ_test_wald <- test_Wald(matrix(univData), p=1, G= 300, Phi = matrix(a_univ), eps=matrix(eps_univ), alpha = .1, estim="DiagH")
univ_test_wald
