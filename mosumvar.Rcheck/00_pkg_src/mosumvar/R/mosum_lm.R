# mosum_lm

## internals
get_H_fixed <- function(mod, X){
  return( -mod$residuals * X)
}

get_A_fixed <- function(H, k,G){
  return(colSums(H[(k+1):(k+G),])-colSums(H[(k-G+1):k,]))
}

get_Tk_fixed <- function(mod, k, G, H, C12 ){ 
  A <- (get_A_fixed(H,k,G)) #difference matrix
  sig <- get_local2(mod,k,G)
  ##Sigma estimator ------
  sgd <- sig^(-.5) * C12
  #------------------------------
  Tkn <-  norm(sgd %*% A,type="2") /sqrt(2*G)#(2*G)^(-.5) 
  return(Tkn)
}

get_Wk_fixed <- function(x, k, G, d, C_root ){ 
  mod1 <- lm(x[(k-G+1):k,1] ~ x[(k-G+1):k,2:d] -1 ) ##left model
  mod2 <- lm(x[(k+1):(k+G),1] ~ x[(k+1):(k+G),2:d] - 1) ##right model
  
  mod1$coefficients[is.na(mod1$coefficients)] <- 0 #handle NAs
  mod2$coefficients[is.na(mod2$coefficients)] <- 0
  
  
  A <- mod2$coefficients - mod1$coefficients
  
  sig <- get_local1(mod1,mod2,G)
  ##Sigma estimator ------
  sgd <- sig^(-.5) * C_root  ##not sqrt of C
  #------------------------------
  Wkn <- sqrt(G/2) *  norm(sgd %*% A,type = "2")#as.numeric( sqrt(t(A) %*% sgd %*% A) )#
  return(Wkn)
}

##variance
get_local1 <- function(mod1,mod2, G){
  out <- (sum(mod1$residuals^2) + sum(mod2$residuals^2))/(2*G) ##local1 variance estim
  return(out)
}

get_local2 <- function(mod,k,G){
  eps <- mod$residuals
  upper <- var(eps[(k+1):(k+G)])
  lower <- var(eps[(k-G+1):(k)])
  out <- (upper + lower)* (G-1)/(2*G)
  return(out)
}

get_global <- function(mod) var(mod$residuals)



## main

#' MOSUM procedure for multivariate regression
#'
#' @param X data matrix with response in column 1, and intercept in any other column
#' @param G integer MOSUM bandwidth
#' @param method string indicating which of `Wald` or `Score` to use
#' @param alpha Numeric significance level
#' @param criterion string location procedure
#' @param nu Numeric location procedure hyperparameter
#' @return list containing Boolean test outcome `Reject`, Numeric rejection threshold `Threshold`, 
#'  Numeric vector of test statistic `mosum`, Integer vector of estimated change points `cps`, Plot `plot`, 
#' @examples
#' data(X0df)
#' mosum_lm(X0df, 200)
#' data(X1df)
#' mosum_lm(X1df, 200)
mosum_lm <- function(X, G, method = c("Wald", "Score")[1], alpha = 0.05, criterion= c("eps","eta")[1], nu=.25){
  X <- as.matrix(X)
  out <- NULL
  n <- dim(X)[1]
  d <- dim(X)[2] 
  #nu <- 0.25
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (d-1)/2 * log(log(n/G)) - log(2/3 * gamma((d-1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector
  C <-t(X[,2:d]) %*% (X[,2:d]) / n
  ev <- eigen(C)
  ##t(ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% (ev$vectors)
  
  if(method == "Wald") { 
    Croot <-  (ev$vectors) %*% diag( (ev$values)^(.5) ) %*% t(ev$vectors) 
    for (t in (G+1):(n-G)) {
      stat[t] <- get_Wk_fixed(X, k=t, G, d = dim(X)[2], Croot ) 
    }
  } 
  if(method == "Score"){  
    C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
    mod <- lm(X[,1] ~ X[,2:d]-1)
    H <- get_H_fixed(mod,X[,2:d])
    for (t in (G+1):(n-G)) {
      stat[t] <- get_Tk_fixed(mod,k=t,G,H,C12)
    }
  }
  
  cps <- get_cps(stat,D_n,G, nu=nu, criterion) #locate
  if(length(cps)>0) Reject<-TRUE
  
  ##Plot------------------------------------
  plot.ts(stat, ylab="Statistic") # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = stat, cps = cps, plot = pl)
  
  return(out)
}



## subsample 
#' MOSUM subsampling procedure for multivariate regression
#'
#' @param X data matrix with response in column 1, and intercept in any other column
#' @param G integer MOSUM bandwidth
#' @param method string indicating which of `Wald` or `Score` to use
#' @param kap Numeric sampling resolution constant
#' @param alpha Numeric significance level
#' @param criterion string location procedure
#' @param nu Numeric location procedure hyperparameter
#' @return list containing Boolean test outcome `Reject`, Numeric rejection threshold `Threshold`, 
#'  Numeric vector of test statistic `mosum`, Integer vector of estimated change points `cps`, Plot `plot`, 
#' @examples
#' data(X0df)
#' mosum_lm_sub(X0df, 200, kap = 1)
#' data(X1df)
#' mosum_lm_sub(X1df, 200, kap = 1)
mosum_lm_sub <- function(X, G, method = c("Wald", "Score")[1], kap = 0.1,  alpha = 0.05, 
                            criterion= c("eps","eta")[1], nu=.25){
  X <- as.matrix(X)
  n <- dim(X)[1]
  d <- dim(X)[2] 
  #nu <- 0.25
  R <- floor(kap*G)
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (d-1)/2 * log(log(n/G)) - log(2/3 * gamma((d-1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  ind <- index(n,G,p=0,kap) #fit grid
  if(n - ind[length(ind)] < G) ind <- head(ind, -1)
  stat <- rep(0, n) #initialise statistic vector
  #if(method == "Wald")statW <- vapply(ind, get_Wkn_RCPP, FUN.VALUE = double(1), x=x,p=p,G=G, estim=estim)#get_W_RCPP(x[ind,],p,G,estim)
  C <-t(X[,2:d]) %*% (X[,2:d]) / n
  ev <- eigen(C)
  ##t(ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% (ev$vectors)
  ##Score Subsample-------------------------------------
  if(method == "Score" || method == "Wald" ){
    stat <- rep(0, n)  
    if(method == "Score"){  
      C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
      for (k in ind) {
        window <- (k-G+1):(k+G)
        mod <- lm(X[window,1] ~ X[window,2:d]-1)
        H <- get_H_fixed(mod,X[window,2:d])
        a <- get_Tk_fixed(mod,k=G,G,H,C12) ##calculate statistic on subsample
        stat[k] <- a ##collect into output vector
        if(a > D_n ){ ##if passes threshold locally
          Reject <- TRUE
          if(k>= 2*G & k <= n-1*G){ #not too close to ends
            #ss <- max(G+1, k-2*G +1); ee <- min(n-G,k+2*G) ##bounds
            newresids <- X[,1] -  as.vector( X[,2:d] %*% mod$coefficients )#predict(mod, as.data.frame(X[ss:ee,2:d]), se.fit=F) #rep(0,  ee-ss) #obtain extended residuals
            mod$residuals <- newresids #overwrite residuals
            newH <-  -1 * newresids *X[,2:d] ##new estimating function series
            max_fill <- min(n-k, 2*G) #for filling in ends
            Tt <- rep(0,  max_fill) #new statistic evaluation
            for (t in  1:max_fill ){
              Tt[t] <- get_Tk_fixed(mod, k =window[t],G,newH, C12)
            }
            
            stat[window[1:max_fill]] <- pmax(stat[window[1:max_fill]], Tt) ##select max at each value
          }
        }  
        
      }
    }
    if(method == "Wald"){  
      statW <- rep(0, length(ind))
      Croot <-  (ev$vectors) %*% diag( (ev$values)^(.5) ) %*% t(ev$vectors) ##matrix root
      for (k in ind) {
        a <- get_Wk_fixed(X,k,G,d,Croot) ##calculate statistic on subsample
        statW[which(ind==k)] <-a ##collect into output vector
      }
      stat[ind] <- statW
      test_stat <- max(stat)
      cps <- c() #empty changepoint vector
      if(test_stat > D_n){ #compare test stat with threshold
        Reject <- TRUE
        C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors) ##need C inverse root for Tk
        times <- which(stat> D_n)#times <- sort(Reduce(union,tlist))
        tlist <- split(times, cumsum(c(1, diff(times) != R))) #split into list of consecutive regions
        
        for (i in 1:length(tlist)) {
          interval <- (max(min(tlist[[i]]) - 2*G, 1)):( min(max(tlist[[i]]) + 2*G, n)) ##extend by +-2G
          #fit var model
          mod <- lm(X[interval,1] ~ X[interval,2:d]-1)
          H <- get_H_fixed(mod,X[interval,2:d])
          
          for (t in (G):(length(interval)-G )) {
            stat[min(interval) -1 + t] <- get_Tk_fixed(mod, t, G, H, C12)
          }
          #stat[interval[(G):(length(interval)-G )]] <- 
          #  vapply( (G):(length(interval)-G ), get_Tk_fixed, FUN.VALUE = 0, mod=mod, G=G,H=H,C12=C12 )  #overwrite statistic
        }
      }
    } 
    
    
    cps <- c()
    if(criterion == "eps"){
      sub_pairs <- get_sub_pairs(stat,D_n,G,kap=kap,nu=nu) #get_sub_pairs
      q <- dim(sub_pairs)[1]
      if(q==0) Reject <- FALSE
      else if (q>0){ ## locate cps
        for (ii in 1:q) {
          interval <- sub_pairs[ii,1]:sub_pairs[ii,2]
          kk <- which.max(stat[interval]) #internal cp location
          cps[ii] <- kk + sub_pairs[ii,1] #- G-p
        }
      }
    }
    if(criterion == "eta"  ){
      cps <- get_local_maxima(stat,D_n,G,nu=2*nu)
      q <- length(cps)
      if(q==0) Reject <- FALSE
    }
  }
  ##BinSeg----------------------------------------
  if(method=="BinSeg"){
    C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
    cps <- MOSUMBS_fixed(X,1,n,D=D_n,G,d,C12,cps=c(), iter=0)
    if(length(cps)==0){ Reject <- F }else{ Reject <- T}
    stat <- c()
  }
  
  ##Plot------------------------------------
  if(method=="Score"||method=="Wald"){
    plot.ts(stat, ylab="Statistic") # plot test statistic
    abline(h = D_n, col = "blue") #add threshold
    if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  }
  if(method=="BinSeg"){
    plot.ts(X[,1], ylab="Y")
    if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  }
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = stat, cps = cps, plot = pl, estim=0)
  
  return(out)
}



## binary segmentation

#' MOSUM Binary Segmentation procedure for multivariate regression
#'
#' @param X data matrix with response in column 1, and intercept in any other column
#' @param G integer MOSUM bandwidth
#' @param alpha Numeric significance level
#' @param max_iter Integer maximum number of Binary Segmentation recursions to allow
#' @return list containing Boolean test outcome `Reject`, Numeric rejection threshold `Threshold`, 
#'  Numeric vector of test statistic `mosum`, Integer vector of estimated change points `cps`, Plot `plot`, 
#' @examples
#' data(X1df)
#' mosum_lm_bs(X1df, 200)
mosum_lm_bs <- function(X, G, alpha = 0.05, max_iter = 10){
  X <- as.matrix(X)
  n <- dim(X)[1]
  d <- dim(X)[2] 
  
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (d-1)/2 * log(log(n/G)) - log(2/3 * gamma((d-1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector
  C <-t(X[,2:d]) %*% (X[,2:d]) / n
  ev <- eigen(C)

  ##BinSeg----------------------------------------
  C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
  cps <- MOSUMBS_fixed(X,1,n,D=D_n,G,d,C12,cps=c(), iter=0)
  if(length(cps)==0){ Reject <- F }else{ Reject <- T}
  stat <- c()
  
  ##Plot------------------------------------
  plot.ts(X[,1], ylab="Y")
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = stat, cps = cps, plot = pl, estim=0)
  
  return(out)
}

MOSUMBS_fixed <- function(x, s, e, D, G, d, C12, cps, iter=0, max_iter = 10){
  if(e-s > 2*G & iter < max_iter){
    iter <- iter + 1 #control recursion depth
    mod <- lm(x[s:e,1] ~ x[s:e,2:d]-1)
    ##eigen
    H <- get_H_fixed(mod,x[s:e,2:d])
    stat <-  vapply((1*G):(e-s-1*G +1), get_Tk_fixed, FUN.VALUE = 0, mod=mod, G=G,H=H,C12=C12 ) 
    TT <- max(stat)
    if(TT > D){ ##passes test
      k <- which.max(stat)+G +s -1 #shift
      cps <- append(cps, k)
      cps1 <- MOSUMBS_fixed(x, s, k, D, G, d, C12, cps, iter)
      cps2 <- MOSUMBS_fixed(x, k+1, e, D, G, d, C12, cps, iter)
      cps <- union(cps, c(cps1,cps2) )
    }
  }
  return( sort(cps) ) #return cps from current iter
}
