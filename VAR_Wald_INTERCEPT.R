## WALD STATISTIC FOR VAR(p) MODEL



##EXAMPLES-------------------------



########################################
##WALD-------------------------------

## Wald H

##estimating fn for channel i, time k
getH_ik_Wald <- function(x, i,k,p, a) { 
  #if(is.null(dim(x))) {x <- matrix(x); Phi <- matrix(Phi); eps <- matrix(eps)} #handle univariate case
  d <- dim(x)[2]
  ai <- a[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1)] #select params for channel i #accounts for intercept
  if(p==1)X <- c(1, x[(k-1),]) #with intercept X <- x[(k-1),]
  if(p==2)X <- c(1, x[(k-1),], x[(k-2),])
  if(p==3)X <- c(1, x[(k-1),], x[(k-2),],x[(k-3),])
  if(p==4)X <- c(1, x[(k-1),], x[(k-2),],x[(k-3),],x[(k-4),])
  y <- as.numeric(x[k,i])
  e <- as.numeric(x[k,i] - ai %*% X) #eps[k,i] #p=1 residual
  H_ik <- t(- y*X + ai%*%X%*%t(X) - e*X) #*2
  return( H_ik)
  #return(list(X%*%t(X), cov(X,X)))
}
#eps_ex <-
#getH_ik_Wald(bf_ts_0, i=1, k=10,p=1, a = (make_a_lu(bf_ts_0, p=1, l= 10, u= 100)))
#getH_ik_Wald(p2_change, i=1, k=10,p=2, a = (make_a_lu(p2_change, p=2, l= 10, u= 100)))

makeH_k_Wald <- function(x, k, p, a){ 
  d <- dim(x)[2]
  if(is.null(dim(x))) {d <- 1} #handle univariate case
  H <- t(rep(0, d+  d*d*p)) #accounts for intercept 
  for (i in 1:d) {
    H[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1)] <- getH_ik_Wald(x,i,k,p,a) #accounts for intercept  H[((i-1)*d+1):( (i-1)*d+d)] <- getH_ik_Wald(x,i,k,a)
  }
  return(t(H))
}
#makeH_k_Wald(bf_ts, k=10, p=1, a = make_a_lu(bf_ts_0, p=1, l= 10, u= 100))
#makeH_k_Wald(p2_change, k=10,p=2, a = (make_a_lu(p2_change, p=2, l= 10, u= 100)))

##evaluate estimating function H locally
makeH_l_u <- function(x, p, l, u, a){ 
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  if(is.null(dim(x))) {n <- length(x); d<- 1} #handle univariate case
  #K <- (G+1):(n-G) #time indices to evaluate over
  H_l_u <- matrix(0, nrow = d+ d*d*p, ncol = u-l+1) #matrix of H values #accounts for intercept
  for (t in 1:(u-l+1)) {
    H_l_u[,t] <- makeH_k_Wald(x, k=l+t-1, p, a) 
  }
  return(H_l_u)
}
#makeH_l_u(x=bf_ts_0, p=1, l=10, u =100, a= (make_a_lu(bf_ts_0, p=1, l= 10, u= 100)))#[,1]
#makeH_l_u(p2_change, p=2, l=10, u=100, a = (make_a_lu(p2_change, p=2, l= 10, u= 100)))

##local (l,u) regression parameter for channel i
get_a_lu_i <- function(x, i,p, l, u){
  y <-   x[l:u,i]  #- mean(x[l:u,i]) #response #intercept
  if(p==1)  {
    Xt <- t(cbind(1, x[(l-1):(u-1),])) #regressors ##p = 1 case only
    y_soln <- colSums(diag(y) %*% t(Xt))
  }  
  if(p==2){
    Xt <- rbind(1, gdata::interleave(t(x[(l-1):(u-1),]),t(x[(l-2):(u-2),]) )) 
    y_soln <- colSums(diag(y) %*% t(Xt))
  }
  if(p==3){
    Xt <- rbind(1, gdata::interleave(t(x[(l-1):(u-1),]),t(x[(l-2):(u-2),]), t(x[(l-3):(u-3),]))) 
    y_soln <- colSums(diag(y) %*% t(Xt))
  }
  if(p==4){
    Xt <- rbind(1, gdata::interleave(t(x[(l-1):(u-1),]),t(x[(l-2):(u-2),]), t(x[(l-3):(u-3),]),t(x[(l-4):(u-4),]) )) 
    y_soln <- colSums(diag(y) %*% t(Xt))
  }  
  #Xt_centred <- Xt - rowMeans(Xt) 
   #t(Xt_centred)
  X_soln <- solve(Xt %*% t(Xt))#cov(t(Xt))
  a_lu_i <- t(y_soln %*% X_soln) 
  return(a_lu_i)
} 
#get_a_lu_i(bf_ts_0, i=1,p=1, 10, 100)
#get_a_lu_i(p2_change, i=1,p=2, 10, 100)

##local (l,u) regression parameter for all channels
make_a_lu <- function(x, p, l, u){
  d <- dim(x)[2]
  a_lu <- rep(0, d + d*d*p) #accounts for intercept
  for (i in 1:d) {
    a_lu[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1)] <- get_a_lu_i(x,i,p,l,u)  #((i-1)*d+1):( (i-1)*d+d)  #accounts for intercept
  }
  return(a_lu)
}
#make_a_lu(bf_ts, p=1, l= 10, u= 100)
#make_a_lu(p2_change,p=2, 10, 100)

## sqrt of local uncentred covariance C_nk
get_V_nk <- function(x, p, l, u){
  d <- dim(x)[2]
  if(p==1)xk <- cbind(1,x[(l-1):(u-1),]) #time sample #intercept
  if(p==2)xk <- cbind(1,x[(l-1):(u-1),],x[(l-2):(u-2),]) 
  if(p==3)xk <- cbind(1,x[(l-1):(u-1),],x[(l-2):(u-2),], x[(l-3):(u-3),], ) 
  if(p==4)xk <- cbind(1,x[(l-1):(u-1),],x[(l-2):(u-2),], x[(l-3):(u-3),], x[(l-4):(u-4),] ) 
  #C <- 1/(u-l) * t(xk) %*% (xk) #uncentred covariance of sample
  #if(d==1) xk <- matrix(xk)
  C <- 1/(u-l) * t(xk) %*% (xk)
  Vlist <- list(1:d)
  for(i in 1:d) {
    Vlist[[i]] <- (C) #coerce into block diagonal form
  }
  out <- Matrix::bdiag(Vlist)
  return(out)
}
#get_V_nk(as.matrix(bf_ts_0), p=1, l= 10, u= 100)
#get_V_nk(univData_0, p=1, l= 10, u= 100)
#get_V_nk(p2_change, p=2, l= 10, u= 100)

## S estimators ----------------------------------------------

## LOCAL1 estimator of variance in channel i at k
getsigma_i_kLOCAL1 <- function(x, i, k, G, a_upper, a_lower) {
  x_upper <- rbind(1,t(x[(k):(k+G-1),])) #upper sample
  res_upper <-  x[(k+1):(k+G), i] - t(a_upper) %*% x_upper #upper residuals
  x_lower <- rbind(1,t(x[(k-G):(k-1),])) #lower sample
  res_lower <-  x[(k-G+1):(k), i] - t(a_lower) %*% x_lower #lower residuals
  sigma_i <- 1/(2*G) * (sum(res_upper^2) + sum(res_lower^2) ) #LOCAL1
  return(sigma_i)
}
#a_upper_example <- get_a_lu_i(bf_ts_0, i=1, l= 100, u= 130)
#a_lower_example <- get_a_lu_i(bf_ts_0, i=1, l= 69, u= 99)
#getsigma_i_kLOCAL1(x = bf_ts_0, i=1, k = 100, G= 30, a_upper = a_upper_example, a_lower = a_lower_example)

getsigma_d_kLOCAL1 <- function(x, k, G,p, a_upper, a_lower){
  d <- dim(x)[2]
  sigma_d <- rep(0, d)
  for (i in 1:d) {
    sigma_d[i] <- getsigma_i_kLOCAL1(x,i,k,G,a_upper[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1)] ,a_lower[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1)]) #accounts for intercept
  }
  return(sigma_d)
}
# a_upper_all_ex <- make_a_lu(bf_ts_0, p=1, l= 100, u= 130)
# a_lower_all_ex <- make_a_lu(bf_ts_0, p=1, l= 69, u= 99)
# bf_sigma_d <- getsigma_d_kLOCAL1(x = bf_ts_0, k = 100, G= 30, a_upper = a_upper_all_ex, a_lower = a_lower_all_ex)

## sigma^2_i estimate for all times k; example function only, NOT used in test
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

get_DiagH_Wald <- function(x, G,p, H_l, H_u){
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  #l <- (k-G):(k-1) #lower time index
  #u <- k:(k+G-1) #upper time index
  H_list <- list(1:d)
  for (i in 1:d) {
    #if(d>1){ #d>1
      H_u_i <- (H_u[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1),] )##LOCALISE ((i-1)*d+1):( (i-1)*d+d) ##accounts for intercept
      H_l_i <- (H_l[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1),] )
      Hbar_u <- rowMeans(H_u_i) 
      H_u_centred <- (H_u_i - Hbar_u)
      Hbar_l <- rowMeans(H_l_i)
      H_l_centred <- (H_l_i - Hbar_l)
      H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred) #sum  } else Hbar_u <- matrix(mean(H_u))
    #} else { #d=1
      #H_u_i <- H_u[,u]
      #H_l_i <- H_l[, l] 
    #  Hbar_u <- mean(H_u)
    #  H_u_centred <- (H_u - Hbar_u) 
    #  Hbar_l <- mean(H_l)
    #  H_l_centred <- (H_l - Hbar_l)
    #  H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred)
    #}
    ##eigen decomposition
    e <- eigen(H_out) #SVD of H_out
    V <- e$vectors
    #if(d>1) {
      H_ <- V %*% diag((e$values)^(-0.5)) %*% t(V)
   # } else {H_ <- e$values^(-0.5) * V %*% t(V)}
    H_list[[i]] <- H_
  }
  Sig_ <-  sqrt(2*G) *Matrix::bdiag(H_list) #coerce into block diagonal form
  return(Sig_)
}
# H_l_ex <- makeH_l_u(x=bf_ts_0, p=1, l=101, u =200, a= (make_a_lu(bf_ts_0, p=1, l= 101, u= 200)))#, eps=bf_ts)
# H_u_ex <- makeH_l_u(x=bf_ts_0, p=1, l=201, u =300, a= (make_a_lu(bf_ts_0, p=1, l= 201, u= 300)))# , eps=bf_ts)
# get_DiagH_Wald(x=bf_ts_0, G=100, H_l = H_l_ex, H_u = H_u_ex)
# H_l_univ <- makeH_l_u(x=univData_0, p=1, l=101, u =200, a= (make_a_lu(univData_0, p=1, l= 101, u= 200)))#, eps=bf_ts)
# H_u_univ <- makeH_l_u(x=univData_0, p=1, l=201, u =300, a= (make_a_lu(univData_0, p=1, l= 201, u= 300)))# , eps=bf_ts)
# get_DiagH_Wald(x=univData_0, G=100, H_l = H_l_univ, H_u = H_u_univ)

get_FullH_Wald <- function(x, G, H_l, H_u){
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  #l <- (k-G):(k-1) #lower time index
  #u <- k:(k+G-1) #upper time index
      Hbar_u <- rowMeans(H_u) 
      H_u_centred <- (H_u - Hbar_u)
      Hbar_l <- rowMeans(H_l)
      H_l_centred <- (H_l - Hbar_l)
      H_out <- (H_l_centred) %*% t(H_l_centred) + (H_u_centred) %*% t(H_u_centred) #sum  } else Hbar_u <- matrix(mean(H_u))

    ##eigen decomposition
    e <- eigen(H_out) #SVD of H_out
    V <- e$vectors
  Sig_ <-  sqrt(2*G) *V %*% diag((e$values)^(-0.5)) %*% t(V) #coerce into block diagonal form
  return(Sig_)
}
#get_FullH_Wald(x=bf_ts_0, G=100, H_l = H_l_ex, H_u = H_u_ex)
#get_FullH_Wald(x=univData_0, G=100, H_l = H_l_univ, H_u = H_u_univ)

get_DiagC_Wald <- function(x, p, sigma_d, k, G){ #Root, NOT Root inverse
  if(is.null(dim(x))) {x <- matrix(x)} #handle univariate case
  n <- dim(x)[1]
  d <- dim(x)[2]
  if(p==1)xk <-  cbind(1, x[(k-G):(k+G-1),]) #lower time sample #includes intercept
  if(p==2)xk <-  cbind(1, x[(k-G):(k+G-1),],x[(k-G-1):(k+G-2),])
  if(p==3)xk <-  cbind(1, x[(k-G):(k+G-1),],x[(k-G-1):(k+G-2),], x[(k-G-2):(k+G-3),] )
  if(p==4)xk <-  cbind(1, x[(k-G):(k+G-1),],x[(k-G-1):(k+G-2),], x[(k-G-2):(k+G-3),], x[(k-G-3):(k+G-4),] )
  #if(d==1) xk <- matrix(xk) #redundant
  C <- 1/(2*G) * t(xk) %*% (xk)
  ##eigen decomposition
  e <- eigen(C) #SVD of C
  V <- e$vectors
  #if(d>1) {
  C_ <- V %*% diag( e$values^(0.5) ) %*% t(V)
  #} else {C_ <- e$values^(-0.5) * V %*% t(V)}
  C_list <- list(1:d)
  for (i in 1:d) {
    C_list[[i]] <- sigma_d[i]^(-0.5) * C_#
  }
  Sig_ <-  (2)^(-0.5) *Matrix::bdiag(C_list) #coerce into block diagonal form
  return(Sig_)
}
#Sig_example <- get_DiagC_Wald(bf_ts_0, bf_sigma_d, 100,30)
#get_DiagC_Wald(matrix(univData), getsigma_d_kLOCAL1(univData,100,100,), 400, 100)
#get_DiagC_Wald(p2_change, getsigma_d_kLOCAL1(p2_change, k=101, G=100,p=2,) , 400, 100)

## W Statistic --------------------------------

## evaluate Wkn (test statistic) for time k
get_Wkn <- function(x, p, k, G, estim){ 
  n<- dim(x)[1]
  d <- dim(x)[2]
  #if(d>1){
  #x_lower_mean <- (colMeans(x[(k-G+1):k,])) ##DYNAMIC CENTERING 
  #x_upper_mean <- (colMeans(x[(k+1):(k+G),]))
  #} else { #d=1
  #  x_lower_mean <- mean(x[(k-G+1):k,]) ##DYNAMIC CENTERING 
  #  x_upper_mean <- mean(x[(k+1):(k+G),])
  #}
  #x_l <- sweep(x,2, x_lower_mean)
  #x_u <- sweep(x, 2, x_upper_mean) 
  a_upper <- make_a_lu(x, p, l=k+1, u=k+G)#x_u
  #res_u <-  x_u - t(a_upper) %*% rbind(rep(0,p),x_u[(p+1):n,]) #upper residuals
  a_lower <- make_a_lu(x, p, l=k-G+1, u=k)#x_l
  #res_l <-  x_l - t(a_lower) %*% rbind(rep(0,p),x_l[(p+1):n,] ) #lower residuals
  ##Sigma estimator options------
  if(estim == "DiagC"){
    sigma_d <- getsigma_d_kLOCAL1(x,k,G,p,a_upper,a_lower)
    Sig_ <- get_DiagC_Wald(x,p,sigma_d,k,G) 
    W_mat <- as.matrix(Sig_ %*% (a_upper-a_lower))
  } else{
    V <- get_V_nk(x, p, k-G+1, k)#x_l
    H_l <- makeH_l_u(x, p, l=k-G+1, u=k , a=a_lower)#x_l
    H_u <- makeH_l_u(x, p, l=k+1, u=k+G , a=a_upper)#x_u
    if(estim == "DiagH")Sig_ <- get_DiagH_Wald(x,G,p,H_l,H_u) #DiagH estimator for Sigma 
    if(estim == "FullH")Sig_ <- get_FullH_Wald(x,G,H_l,H_u) #FullH estimator 
    W_mat <- as.matrix(Sig_ %*%V %*% (a_upper-a_lower)) #argument for norm
  }
  
  #------------------------------
  #W <- W + 1/sigma_i * norm(root_C %*% (a_i_upper - a_i_lower), type = "F")
 
  W <- sqrt(G/2) * norm(W_mat, type="F")
  return(W)
}
#get_Wkn(bf_ts_0, p=1, k= 101, G=100, estim = "FullH")
#get_Wkn(bf_ts_0, p=1, k= 101, G=100, estim = "DiagC")
#get_Wkn(univData_0, p=1, k= 101, G=100, estim = "DiagH")
#get_Wkn(p2_change, p=2, k= 105, G=100, estim = "DiagH")

get_W <- function(x, p, G, estim){
  n <- dim(x)[1]
  K <- (G+p):(n-G)
  #H_all <- makeH_all(x, p, G, Phi, eps) ##WRONG - use local H
  out <- rep(0, n)
  for (k in K) {
    out[k] <- get_Wkn(x,p,k,G,estim)
  }
  return(out)
}
#W_ex <- get_W(x= bf_ts_0, p=1, G=100, estim="DiagH")
#get_W(x= p2_change, p=2, G=100, estim="DiagH")

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
change_test_wald_C <- test_Wald(x= var_change_0, p=1, G= 200, alpha = .1, estim = "DiagC")
change_test_wald_C$plot
# no change
nochange_test_wald <- test_Wald(x=nochange_0, p=1, G=150, alpha = 0.1) 
nochange_test_wald$plot
nochange_test_wald_C <- test_Wald(x=nochange_0, p=1, G=150, alpha = 0.1,estim = "DiagC") 
#univariate
univ_test_wald <- test_Wald(univData_0, p=1, G= 300, alpha = .1, estim="DiagH")
univ_test_wald$plot
univ_test_wald_C <- test_Wald(univData_0, p=1, G= 300, alpha = .1, estim="DiagC")

# p=2 change
p2_test_wald <- test_Wald(p2_change, p=2, G= 300, alpha = .1, estim="DiagH")

## benchmark---------------------------
library(microbenchmark)
library(ggplot2)
dcw <- function()test_Wald(x=var_change, p=1, G=200, alpha = 0.05, estim = "DiagC") 
dhw <- function()test_Wald(x=var_change, p=1, G=200, alpha = 0.05, estim = "DiagH") 
fhw <- function()test_Wald(x=var_change, p=1, G=200, alpha = 0.05, estim = "FullH") 
        
mbw <- microbenchmark(dcw(),dhw(),fhw(), times = 10)
autoplot(mbw) + ggtitle("Wald Procedure Runtime, by Estimator") 
