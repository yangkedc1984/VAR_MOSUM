## SCORE STATISTIC FOR VAR(p) MODEL


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
##SCORE----------------------------------

##estimating fn for channel i, time k
getH_ik <- function(x, i,k, Phi, eps) { 
  a <- t(Phi[i,]) #A_1[i,]
  X <- x[(k-1),]
  #X <- t(X)
  y <- as.numeric(x[k,i])
  e <- eps[k,i]
  H_ik <- - y*X + a%*%X%*%t(X) - e*X
  return(t(H_ik))
  #return(list(X%*%t(X), cov(X,X)))
}
#getH_ik(bf_ts, i=1, k=10, Phi = A_1, eps = eps)

##make whole estimating function H (pd^2 vector) at time k
makeH_k <- function(x, k, p, Phi, eps){ 
  d <- dim(x)[2]
  H <- t(rep(0, d*d*p))
  for (i in 1:d) {
    H[((i-1)*d+1):( (i-1)*d+d)] <- getH_ik(x,i,k,Phi,eps)
  }
  return(t(H))
}
#makeH_k(bf_ts, k=10, p=1, Phi = A_1, eps = eps)

##evaluate estimating function H at all time steps
makeH_all <- function(x, p, G, Phi, eps){ 
  #wrap_makeH_k <- function(k) makeH_k(x, k, p, Phi, eps) #wrapper function for makeH_k
  n <- dim(x)[1]
  d <- dim(x)[2]
  K <- (G+1):(n-G) #time indices to evaluate over
  H_all <- matrix(0, nrow = d*d*p, ncol = n) #matrix of H values
  for (t in K) {
    H_all[,t] <- makeH_k(x, k=t, p, Phi, eps) 
  }
  return(H_all)
}
#H_all_example <- makeH_all(x=bf_ts, p=1, G=5, Phi=A_1, eps=eps)

##Difference Vector at k
getA <- function(x, k, p, G, Phi, eps, H_all){ 
  r <- H_all[,(k+1):(k+G)]; l <- H_all[,(k-G+1):(k)] #left/right of window
  A <- rowSums(r) - rowSums(l)#difference
  A <- as.matrix(A)
  return(A)
}
#getA(x=bf_ts, k = 10, p=1, G=5, Phi=A_1, eps=eps, H_all = H_all_example)


##sigma_i options---------------------------------------------
##global estimate for sigma^2_i 
getsigma_iGlobal <- function(eps,i){
  mean(eps[-1,i]^2)
}
#getsigma_iGlobal(eps,1)

#getsigma_iLOCAL1 <- function(eps,i)

##estimate inverse sqrt of long-run covariance, Sigma^{-1/2}
getSigma_root_inv <- function(x, eps, p){ 
  n <- dim(x)[1]
  d <- dim(x)[2]
  sigma_d <- rep(0,d)
  for (i in 1:d) { #sigma(i)^{-1/2}
    sigma_d[i] <- getsigma_iGlobal(eps,i)^(-0.5)
  }
  
  c <- 1/n * t(x) %*% x #c <- cov(x)
  e <- eigen(c) #SVD of c
  V <- e$vectors
  c_ <- V %*% diag((e$values)^(-.5)) %*% t(V) #c^{1-/2}
  c_list <- list(1:d)
  for (i in 1:d) {
    c_list[[i]] <- sigma_d[i]*c_
  }
  Sig_ <-  Matrix::bdiag(c_list) #coerce into bloack diagonal form
  return(Sig_)
  #return(c_list)
}
#Sig_example <- getSigma_root_inv(bf_ts,eps,1)

##evaluate statistic at time k
getTkn <- function(x, k, p, G, Phi, eps, H_all, Sig_){ 
  A <- getA(x, k, p, G, Phi, eps, H_all) #difference matrix
  
  T_in <- as.matrix(Sig_%*%A)
  Tkn <- (2*G)^(-.5) * norm(T_in,type = "2")
  return(Tkn)
}
#getTkn(x=bf_ts, k=100, p=1, G=20, Phi=A_1, eps=eps, H_all = H_all_example, Sig_ = Sig_example)

##get statistic over all times K
getT <- function(x, p, G, Phi, eps){ 
  n <- dim(x)[1]
  K <- (G+1):(n-G) #time indices to evaluate over
  H_all <- makeH_all(x, p, G, Phi, eps)
  Tkn <- rep(0,n)
  Sig_ <- getSigma_root_inv(x,eps,p) #Inverse root long-run covariance
  for(k in K){
    Tkn[k] <- getTkn(x,k,p,G,Phi,eps, H_all, Sig_) 
  }
  #return(Tkn)
  return(Sig_)
}
#Tn_example <- getT(x=bf_ts, p=1, G= 20, Phi=A_1, eps=eps)
#getT(x=var_change, p=1, G=40, Phi = a_change, eps = eps_change)

##get change point estimates
get_cps <- function(Tn, D_n, nu = 1/4){
    n <- length(Tn)
    lshift <- c(Tn[-1],0); rshift <- c(0,Tn[-n]); 
    over <- (Tn >D_n) #indices are greater than D_n?
    v <- which(over && (lshift < D_n) )#lowers
    v <- c(1,v) #append n
    w <- which(over && (rshift < D_n) )#uppers
    w <- c(w,n) #append 0
    q <- length(w) #number of CPs
    cps <- rep(0, q)
    for (i in 1:q) {
      cps[i] <- v[i] + which.max(Tn[ (v[i]):(w[i]) ] )
    }
    return(cps)
}
#get_cps(Tn = Tn_example, D_n = 5,nu = 1/4)

##Score-type test
test_Score <- function(x, p, G, Phi, eps, alpha = 0.05){ 
  n <- dim(x)[1] #dimensions
  d <- dim(x)[2] 
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d/2 * log(log(n/G)) - log(2/3 * gamma(d/2))
  D_n <- (b+c_alpha)/a #threshold
  Reject <- FALSE
  ##Run test-----------------------------
  Tn <- ts(getT(x,p,G,Phi,eps)) #evaluate statistic at each time k
  test_stat <- max(Tn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Tn,D_n, nu=1/4)
  } 
  ##Plot------------------------------------
  plot(Tn) # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = Tn, cps = cps)
  return(out)
}


##  bf example
test_Score(x=bf_ts, p=1, G=150, Phi = A_1, eps = eps, alpha = 0.1)
##  change example
test_Score(x=var_change, p=1, G=150, Phi = a_change, eps = eps_change, alpha = 0.1) 




