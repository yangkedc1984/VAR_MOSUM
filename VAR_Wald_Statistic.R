## WALD STATISTIC FOR VAR(p) MODEL



##EXAMPLES-------------------------
##example------------------------- bankfox data
bf_ts <- ts(bf[672:dim(bf)[1], c(1, 3:6, 8:10)])
bf_ts <- na.fill(bf_ts, fill = "extend")
ts.plot(bf_ts)

model <- ar(bf_ts, demean = F, method = "ols")
#Phi <- as.vector(t(as.matrix(model$ar)))
A_1 <- matrix(model$ar, nrow=8, ncol=8)
eps <- model$resid

##example------------------------- change data
rData1 <- rSim(a1, e1)
rData1[1,] <- runif(5, 0, 0.02) #prevent NaN
rData2 <- rSim(a2, e2)
rData2[1,] <- runif(5, 0, 0.02)


var_change <- ts(rbind(rData1, rData2))
change_model <- ar(var_change, order.max = 1, demean = F, method = "ols")
a_change <- matrix(change_model$ar, nrow=5, ncol=5)
eps_change <- change_model$resid
plot(var_change)


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
get_root_C_nk <- function(x, p, l, u){
  n <- dim(x)[2]
  xk <- x[(l-1):(u-1),] #time sample
  C <- 1/n * t(xk) %*% (xk) #uncentred covariance of sample
  e <- eigen(C) #SVD of C
  V <- e$vectors
  C_ <- V %*% diag((e$values)^(.5)) %*% t(V) #c^{1/2}
  return(C_)
}
#get_root_C_nk(bf_ts, p=1, l= 10, u= 100)

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

#get_Gamma <- function(){
#}

get_Wkn <- function(x, p, k, G){
  d <- dim(x)[2]
  root_C <- get_root_C_nk(x, p, k-G+1, k+G)
  W <- 0 #initialise stat
  for (i in 1:d) {
    a_i_upper <- get_a_lu_i(x, i, l=k, u=k+G)
    a_i_lower <- get_a_lu_i(x, i, l=k-G-1, u=k-1)
    sigma_i <- sqrt( getsigma_i_kLOCAL1(x, i, k, G, a_upper = a_i_upper, a_lower = a_i_lower) )
    W <- W + 1/sigma_i * norm(root_C %*% (a_i_upper - a_i_lower), type = "2")
  }
  W <- sqrt(G) * W
  return(W)
}
get_Wkn(x = bf_ts, p=1, k= 33, G=30)

get_W <- function(x, p, G){
  n <- dim(x)[1]
  K <- (G+3):(n-G-2)
  out <- rep(0, n)
  for (k in K) {
    out[k] <- get_Wkn(x,p,k,G)
  }
  #out <- sapply(K, get_Wkn, x=x, p=p, G=G)
  #zeros <- rep(0,G)
  #cbind(zeros, out, zeros)
  return(out)
}
#get_W(x=bf_ts, p=1, G=30)

##get change point estimates -- same as Score procedure
#get_cps <- function(Tn, D_n, nu = 1/4){
#  n <- length(Tn)
#  lshift <- c(Tn[-1],0); rshift <- c(0,Tn[-n]); 
#  over <- (Tn >D_n) #indices are greater than D_n?
#  v <- which(over && (lshift < D_n) )#lowers
#  v <- c(1,v) #append n
#  w <- which(over && (rshift < D_n) )#uppers
#  w <- c(w,n) #append 0
#  q <- length(w) #number of CPs
#  cps <- rep(0, q)
#  for (i in 1:q) {
#    cps[i] <- v[i] + which.max(Tn[ (v[i]):(w[i]) ] )
#  }
#  return(cps)
#}

##Score-type test
test_Wald <- function(x, p, G, alpha = 0.05){ 
  n <- dim(x)[1] #dimensions
  d <- dim(x)[2] 
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d/2 * log(log(n/G)) - log(2/3 * gamma(d/2))
  D_n <- (b+c_alpha)/a #threshold
  Reject <- FALSE
  ##Run test-----------------------------
  Wn <- ts(get_W(x,p,G)) #evaluate statistic at each time k
  test_stat <- max(Wn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Wn,D_n, nu=1/4)
  } 
  ##Plot------------------------------------
  plot(Wn) # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = Wn, cps = cps)
  return(out)
}

## examples
# bf
test_Wald(x=bf_ts, p=1, G=100, alpha= .05)

# change
test_Wald(x= var_change, p=1, G= 300, alpha = .1)
