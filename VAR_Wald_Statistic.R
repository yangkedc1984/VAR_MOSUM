## WALD STATISTIC FOR VAR(p) MODEL



##EXAMPLES-------------------------



########################################
##WALD-------------------------------

## Wald H

##estimating fn for channel i, time k
getH_ik_Wald <- function(x, i,k, a) { 
  #if(is.null(dim(x))) {x <- matrix(x); Phi <- matrix(Phi); eps <- matrix(eps)} #handle univariate case
  d <- dim(x)[2]
  ai <- a[((i-1)*d+1):( (i-1)*d+d)] #select params for channel i
  X <- x[(k-1),]
  #X <- t(X)
  y <- as.numeric(x[k,i])
  e <- as.numeric(x[k,i] - ai %*% X) #eps[k,i] #p=1 residual
  H_ik <- t(- y*X + ai%*%X%*%t(X) - e*X) #*2
  return( H_ik)
  #return(list(X%*%t(X), cov(X,X)))
}
#eps_ex <-
getH_ik_Wald(bf_ts_0, i=1, k=10, a = (make_a_lu(bf_ts_0, p=1, l= 10, u= 100)))

makeH_k_Wald <- function(x, k, p, a){ 
  d <- dim(x)[2]
  if(is.null(dim(x))) {d <- 1} #handle univariate case
  H <- t(rep(0, d*d*p))
  for (i in 1:d) {
    H[((i-1)*d+1):( (i-1)*d+d)] <- getH_ik_Wald(x,i,k,a)
  }
  return(t(H))
}
makeH_k_Wald(bf_ts, k=10, p=1, a = make_a_lu(bf_ts_0, p=1, l= 10, u= 100))

##evaluate estimating function H locally
makeH_l_u <- function(x, p, l, u, a){ 
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  if(is.null(dim(x))) {n <- length(x); d<- 1} #handle univariate case
  #K <- (G+1):(n-G) #time indices to evaluate over
  H_l_u <- matrix(0, nrow = d*d*p, ncol = u-l+1) #matrix of H values
  for (t in 1:(u-l+1)) {
    H_l_u[,t] <- makeH_k_Wald(x, k=l+t-1, p, a) 
  }
  return(H_l_u)
}
makeH_l_u(x=bf_ts_0, p=1, l=10, u =100, a= (make_a_lu(bf_ts_0, p=1, l= 10, u= 100)))#[,1]


##local (l,u) regression parameter for channel i
get_a_lu_i <- function(x, i, l, u){
  y <- x[l:u,i] #- mean(x[l:u,i]) #response
  Xt <- t(x[(l-1):(u-1),]) #regressors ##p = 1 case only
  #Xt_centred <- Xt - rowMeans(Xt) 
  y_soln <- colSums(diag(y) %*% t(Xt)) #t(Xt_centred)
  X_soln <- solve(Xt %*% t(Xt))#cov(t(Xt))
  a_lu_i <- t(y_soln %*% X_soln) 
  return(a_lu_i)
} 
get_a_lu_i(bf_ts_0, 1, 10, 100)

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
  #C <- 1/(u-l) * t(xk) %*% (xk) #uncentred covariance of sample
  if(d==1) xk <- matrix(xk)
  C <- cov(xk)
  Vlist <- list(1:d)
  for(i in 1:d) {
    Vlist[[i]] <- (C) #coerce into block diagonal form
  }
  out <- Matrix::bdiag(Vlist)
  return(out)
}
#get_V_nk(as.matrix(bf_ts_0), p=1, l= 10, u= 100)
#get_V_nk(univData_0, p=1, l= 10, u= 100)

## LOCAL1 estimator of variance in channel i at k
getsigma_i_kLOCAL1 <- function(x, i, k, G, a_upper, a_lower) {
  x_upper <- t(x[(k):(k+G-1),]) #upper sample
  res_upper <-  x[(k+1):(k+G), i] - t(a_upper) %*% x_upper #upper residuals
  x_lower <- t(x[(k-G):(k-1),]) #lower sample
  res_lower <-  x[(k-G+1):(k), i] - t(a_lower) %*% x_lower #lower residuals
  sigma_i <- 1/(2*G) * (mean(res_upper^2) +mean(res_lower^2) ) #LOCAL1
  return(sigma_i)
}
#a_upper_example <- get_a_lu_i(bf_ts_0, i=1, l= 100, u= 130)
#a_lower_example <- get_a_lu_i(bf_ts_0, i=1, l= 69, u= 99)
#getsigma_i_kLOCAL1(x = bf_ts_0, i=1, k = 100, G= 30, a_upper = a_upper_example, a_lower = a_lower_example)

getsigma_d_kLOCAL1 <- function(x, k, G, a_upper, a_lower){
  d <- dim(x)[2]
  sigma_d <- rep(0, d)
  for (i in 1:d) {
    sigma_d[i] <- getsigma_i_kLOCAL1(x,i,k,G,a_upper[((i-1)*d+1):( (i-1)*d+d)] ,a_lower[((i-1)*d+1):( (i-1)*d+d)])
  }
  return(sigma_d)
}
#a_upper_all_ex <- make_a_lu(bf_ts_0, p=1, l= 100, u= 130)
#a_lower_all_ex <- make_a_lu(bf_ts_0, p=1, l= 69, u= 99)
#getsigma_d_kLOCAL1(x = bf_ts_0, k = 100, G= 30, a_upper = a_upper_all_ex, a_lower = a_lower_all_ex)

## sigma^2_i estimate for all times k; example function only, not used in test
getsigma_i_all <- function(x,i,G, a_upper, a_lower){
  n <- dim(x)[1]
  d <- dim(x)[2]
  K <- (G+3):(n-G-2)
  a_i_upper <- a_upper[((i-1)*d+1):( (i-1)*d+d)] #get_a_lu_i(x, i, l=k+1, u=k+G) ##WRONG - should be local
  a_i_lower <- a_lower[((i-1)*d+1):( (i-1)*d+d)]#get_a_lu_i(x, i, l=k-G+1, u=k)
  
  out <- rep(0, n)    
  for (k in K) {
    out[k] <- getsigma_i_kLOCAL1(x,i,k,G, a_i_upper, a_i_lower)
  }
  return(ts(out))
}
#getsigma_i_all(x=bf_ts_0,i=1,G=50)
#plot(1/sqrt(getsigma_i_all(x=bf_ts,i=1,G=50)) )
#plot(1/sqrt(getsigma_i_all(x=var_change,i=1,G=50)) ); abline( h= getsigma_iGlobal(eps,1)^(-0.5) )

## Sigma estimators - Wald Specific -------------------------------

get_DiagH_Wald <- function(x, G, H_l, H_u){
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  #l <- (k-G):(k-1) #lower time index
  #u <- k:(k+G-1) #upper time index
  H_list <- list(1:d)
  for (i in 1:d) {
    if(d>1){ #d>1
      H_u_i <- (H_u[((i-1)*d+1):( (i-1)*d+d),] )##LOCALISE
      H_l_i <- (H_l[((i-1)*d+1):( (i-1)*d+d),] )
      Hbar_u <- rowMeans(H_u_i) 
      H_u_centred <- (H_u_i - Hbar_u)
      Hbar_l <- rowMeans(H_l_i)
      H_l_centred <- (H_l_i - Hbar_l)
      H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred) #sum  } else Hbar_u <- matrix(mean(H_u))
    } else { #d=1
      #H_u_i <- H_u[,u]
      #H_l_i <- H_l[, l] 
      Hbar_u <- mean(H_u)
      H_u_centred <- (H_u - Hbar_u) 
      Hbar_l <- mean(H_l)
      H_l_centred <- (H_l - Hbar_l)
      H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred)
    }
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
H_l_ex <- makeH_l_u(x=bf_ts_0, p=1, l=101, u =200, a= (make_a_lu(bf_ts_0, p=1, l= 101, u= 200)))#, eps=bf_ts)
H_u_ex <- makeH_l_u(x=bf_ts_0, p=1, l=201, u =300, a= (make_a_lu(bf_ts_0, p=1, l= 201, u= 300)))# , eps=bf_ts)
get_DiagH_Wald(x=bf_ts_0, G=100, H_l = H_l_ex, H_u = H_u_ex)
H_l_univ <- makeH_l_u(x=univData_0, p=1, l=101, u =200, a= (make_a_lu(univData_0, p=1, l= 101, u= 200)))#, eps=bf_ts)
H_u_univ <- makeH_l_u(x=univData_0, p=1, l=201, u =300, a= (make_a_lu(univData_0, p=1, l= 201, u= 300)))# , eps=bf_ts)
get_DiagH_Wald(x=univData_0, G=100, H_l = H_l_univ, H_u = H_u_univ)


## W Statistic --------------------------------

## evaluate Wkn (test statistic) for time k
get_Wkn <- function(x, p, k, G, estim){ 
  n<- dim(x)[1]
  d <- dim(x)[2]
  if(d>1){
  x_lower_mean <- (colMeans(x[(k-G+1):k,])) ##DYNAMIC CENTERING 
  x_upper_mean <- (colMeans(x[(k+1):(k+G),]))
  } else { #d=1
    x_lower_mean <- mean(x[(k-G+1):k,]) ##DYNAMIC CENTERING 
    x_upper_mean <- mean(x[(k+1):(k+G),])
  }
  x_l <- sweep(x,2, x_lower_mean)
  x_u <- sweep(x, 2, x_upper_mean) 
  V <- get_V_nk(x_l, p, k-G+1, k)
    
  a_upper <- make_a_lu(x_u, p, l=k+1, u=k+G)
  #res_u <-  x_u - t(a_upper) %*% rbind(rep(0,p),x_u[(p+1):n,]) #upper residuals
  a_lower <- make_a_lu(x_l, p, l=k-G+1, u=k)
  #res_l <-  x_l - t(a_lower) %*% rbind(rep(0,p),x_l[(p+1):n,] ) #lower residuals
  ##Sigma estimator options------
  if(estim == "DiagC"){
    #sigma_d <- getsigma_d_kLOCAL1(x,k,G,a_upper,a_lower)
    #Sig_ <- get_DiagC_rootinv(x,eps,sigma_d,k,G) ## NOT WALD-SPECIFIC
  } else{
    H_l <- makeH_l_u(x_l, p, l=k-G+1, u=k , a=a_lower)
    H_u <- makeH_l_u(x_u, p, l=k+1, u=k+G , a=a_upper)
    if(estim == "DiagH")Sig_ <- get_DiagH_Wald(x,G,H_l,H_u) #DiagH estimator for Sigma 
    #if(estim == "FullH")Sig_ <- get_FullH_rootinv(x,k,G,H_all) #FullH estimator ## NOT WALD-SPECIFIC
  }
  
  #------------------------------
  #W <- W + 1/sigma_i * norm(root_C %*% (a_i_upper - a_i_lower), type = "F")
  W_mat <- as.matrix(Sig_ %*%V %*% (a_upper-a_lower)) #argument for norm
  W <- sqrt(G/2) * norm(W_mat, type="F")
  return(W)
}
get_Wkn(bf_ts_0, p=1, k= 101, G=100, estim = "DiagH")
get_Wkn(univData_0, p=1, k= 101, G=100, estim = "DiagH")

get_W <- function(x, p, G, estim){
  n <- dim(x)[1]
  K <- (G+1):(n-G)
  #H_all <- makeH_all(x, p, G, Phi, eps) ##WRONG - use local H
  out <- rep(0, n)
  for (k in K) {
    out[k] <- get_Wkn(x,p,k,G,estim)
  }
  return(out)
}
W_ex <- get_W(x= bf_ts_0, p=1, G=100, estim="DiagH")

##get change point estimates -- same as Score procedure

## TEST ----------------------------------------------------
##Wald-type test
test_Wald <- function(x, p, G, alpha = 0.05, estim="DiagH"){ 
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
  Wn <- ts(get_W(x,p,G,estim)) #evaluate statistic at each time k
  test_stat <- max(Wn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Wn,D_n,G, nu=1/4)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
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


## examples-------------------------------------------------------------
# bf
bf_test_wald <- test_Wald(x= bf_ts_0, p=1, G=200, alpha= .1, estim = "DiagH")
bf_test_wald_C <- test_Wald(x= bf_ts_0, p=1, G=200, Phi = A_1, eps=eps, alpha= .1, estim = "DiagC")
bf_test_wald$plot
# change
change_test_wald <- test_Wald(x= var_change_0, p=1, G= 200, alpha = .1, estim = "DiagH")
change_test_wald$plot
# no change
nochange_test_wald <- test_Wald(x=nochange_0, p=1, G=150, alpha = 0.1) 
nochange_test_wald$plot
#univariate
univ_test_wald <- test_Wald(univData_0, p=1, G= 300, alpha = .1, estim="DiagH")
univ_test_wald$plot
