# ## Subsampling algorithm
# library(Rcpp)
# library(RcppParallel)
# library(RcppArmadillo)
# library(Matrix)
# sourceCpp(file = "Wald_RcppParallel.cpp")
# sourceCpp(file = "Score_Rcpp.cpp")
# ##

###################
# WARNING MESSAGE #
###################

dim_warning <- function(n, G, d, p, method) {
  dim <- d*(d*p + 1 ) + d*(d+1)/2
  dimScore <- d*(d*p + 1 )/2 + d*(d+1)/2
  fl <- floor(n^(2/3))


  W3dim <- paste0("Bandwidth too small relative to model dimensions: set G > d(dp + 1) * log(d(dp + 1)) = ", d*(d*p + 1)* log(d(dp + 1)),
                  "\n")
  Wfl <- paste0("Bandwidth too small relative to sample size: set G > floor(n^(2/3)) = ", fl, "\n" )
  Wlarge <- "Large dimensions: consider `option = univariate`\n"


  if(G < dim & method == "Wald") warning(paste0("Not enough degrees of freedom for Wald method: set G > d(dp + 1) + d(d+1)/2 = ", dim, "\n"))
  if(G < dimScore & method == "Score")warning(paste0("Not enough degrees of freedom for Score method: set G > d(dp + 1)/2 + d(d+1)/2 = ", dimScore, "\n"))
  if(G < 3*dim) warning(W3dim)
  if(G < fl ) warning(Wfl)
  if(d*(d*p + 1 ) > 30) warning(Wlarge)

}
# dim_warning(100, 10, 5,5)
# dim_warning(100, 145, 5,5)


index <- function(n, G, p, kap){ #define grid
  a <- 1:n
  R <- floor(kap*G)
  b <- a[seq.int(G+p, n-G, R)]
  return(b)
}
#index(1000, 200, 2, 0.1)

##get subsample significant pairs
get_sub_pairs <- function(Tn, D_n, G, kap, nu = 1/4){
  n <- length(Tn)
  R <- 1#floor(kap*G)
  rshift <- c(Tn[-(1:R)],rep(0,R)); lshift <- c(rep(0,R),Tn[- ((n-R+1):n)]);
  over <- (Tn >D_n) #indices greater than D_n
  v <- which(over & lshift < D_n )#& lshift >0 ) #lowers
  w <- which(over & rshift < D_n )#& rshift >0) #uppers
  nu_remove <- which(w-v >= nu*R*G) #nu(epsilon) test for distance between
  v_nu <- v[nu_remove]; w_nu <-w[nu_remove] #c(w[nu_remove],n)
  sub_pairs <- cbind(v_nu,w_nu)
  return(sub_pairs)
  #plot(lshift); lines(Tn)
}
#get_sub_pairs(Tn = msub$mosum, D_n = 5, G= 200, kap = 0.3, nu = .25)

get_local_maxima <- function(Tn, D_n, G, nu = 1/4) {
  n <- length(Tn)
  cps <- c()
  window <- floor(nu*G)
  for(t in (G+1):(n-G)){#all(Tn[(t - window):(t+window)]  > D_n) &
    if( Tn[t] > D_n & Tn[t] == max(Tn[(t - window):(t+window)]) ){cps <- append(cps,t)} ##add to list
  }
  return(cps)
}
#get_local_maxima(Tn = msub$mosum, D_n = 3, G= 200,  nu = .25)

##get subsample change point estimates

mosum_sub <- function(x, p, G, method = "Wald", estim = "DiagC", varEstim = "Local", kap = 1,  alpha = 0.05, criterion="eps", nu=.25){
  n <- dim(x)[1]
  d <- dim(x)[2]

  dim_warning(n,G,d,p, method)

  R <- floor(kap*G)
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d*(d*p+1)/2 * log(log(n/G)) - log(2/3 * gamma(d*(d*p+1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  ind <- index(n,G,p,kap) #fit grid
  stat <- rep(0, n) #initialise statistic vector
  if(method == "Wald"){
    statW <- vapply(ind, get_Wkn_RCPP, FUN.VALUE = double(1), x=x,p=p,G=G, estim=estim)#get_W_RCPP(x[ind,],p,G,estim)
    stat[ind] <- statW
    #return(stat)
    test_stat <- max(stat)
    cps <- c() #empty changepoint vector
    if(test_stat > D_n){ #compare test stat with threshold
      Reject <- TRUE
      times <- which(stat> D_n)#times <- sort(Reduce(union,tlist))
      tlist <- split(times, cumsum(c(1, diff(times) != R))) #split into list of consecutive regions

      for (i in 1:length(tlist)) {
        interval <- (max(min(tlist[[i]]) - 2*G, 1)):( min(max(tlist[[i]]) + 2*G, n)) ##extend by +-2G
        #fit var model
        mod <- ar(x[interval,], order.max = p, demean = T, method = "ols", aic = F)
        mod_a <- mod$x.intercept
        eps <- mod$resid; eps[1:p,] <- 1e-4 ##solve NA
        if(p==1) mod_a <- cbind(mod_a,  matrix( mod$ar, nrow=d, ncol=d))
        if(p>1){
          for (jj in 1:p){ #collect parameters into mat
            mod_a <- cbind(mod_a,  mod$ar[jj,,])
          }
        }
        #interval.sort <- sort( union(interval-p,interval))
        #stat[(min(tlist[[i]])):(max(tlist[[i]]))] <- get_T_RCPP( as.matrix(x[interval,]),p,G,Phi= as.matrix(mod_a), eps=as.matrix(eps),estim = estim) #overwrite statistic
        stat[interval[(1*G):(length(interval)-1*G )]] <- get_T_RCPP(
          as.matrix(x[interval,]),p,G,Phi= as.matrix(mod_a), eps=as.matrix(eps),estim = estim)[(1*G):(length(interval)-1*G )] #overwrite statistic
      }
    }
  }

  if(method == "Score"){
    #statW <- rep(0, length(ind))
    for (k in ind){
      s <- max(k-G, 1); e <- min(n,k+G) #start and end of intervals
      subsample <- s:e #set subsample indices
      mod <- ar.ols(x[subsample,],aic=F,order.max = p,demean = T)
      ##build parameter matrix
      if(p==1)a <- cbind(mod$x.intercept, matrix( mod$ar, nrow=d, ncol=d))
      if(p>1) {
        a <- mod$x.intercept
        for (pp in 1:p) {
          a <- cbind(a,matrix( mod$ar[pp,,], nrow=d, ncol=d) ) ##append lagged parameters
        }
      }
      eps <- mod$resid; eps[1:p,] <- 1e-4 ##solve NA
      Tk <- get_T_RCPP(x[subsample,], p,G-p-1,a,eps,estim,var_estim = varEstim) ##calculate statistic on subsample
      stat[k] <- max(Tk) ##collect into output vector
      if(stat[k] > D_n ){ ##if passes threshold locally
        Reject <- TRUE
        if(k> 2*G & k <= n-1*G){ #not too close to ends
          ss <- max(G+p+1, s-G+1-p); ee <- min(n,e+G) ##bounds
          newresids <- matrix(0, nrow = ee-ss, ncol = d) #obtain extended residuals
          for(t in 1:(ee-ss) ){
            newresids[t,] <- predict(mod,x[ss - 1 + (t-p):t,], se.fit=F)
          }
          Tt <-  get_T_RCPP(x[ss:ee,], p,G-p-1,a,newresids,estim,var_estim = varEstim) ##evaluate G-window locally
          stat[subsample] <- pmax(stat[subsample], Tt[G:(3*G)]) ##select max at each value
        }
      }
    }
  }

  cps1 <- cps2 <- c() ##assign empty cps
  if(criterion %in% c("eps","union") ){
    sub_pairs <- get_sub_pairs(stat,D_n,G,kap=kap,nu=nu) #get_sub_pairs
    q <- dim(sub_pairs)[1]
    if(q==0) Reject <- FALSE
    else if (q>0){ ## locate cps
      for (ii in 1:q) {
        interval <- sub_pairs[ii,1]:sub_pairs[ii,2]
        kk <- which.max(stat[interval]) #internal cp location
        cps1[ii] <- kk + sub_pairs[ii,1] #- G-p
      }
    }
  }
  if(criterion %in% c("eta","union") ){
    cps2 <- get_local_maxima(stat,D_n,G,nu=2*nu)
    q <- length(cps2)
    if(q==0) Reject <- FALSE
  }
  cps <- union(cps1,cps2) ##output union

  ##Plot------------------------------------
  plot.ts(stat, ylab="Statistic") # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = stat, cps = cps, plot = pl, estim=estim)
  return(out)
}

# msub <- mosum_sub(x=p2_change,p=2,G=200, method = "Wald", estim = "DiagC", kap = 0.6)
# mmm <- mosum_sub(x=dp2_change,p=2,G=200, estim = "DiagC", kap = 1)
#   test_Wald_new(x=dp2_change,p=2,G=200, alpha = 0.05,  estim = "DiagC")
#
#
# msubScore <- mosum_sub(x=dp2_change,p=2,G=200, method = "Score", estim = "DiagC", kap = 0.6)
#
#
# sourceCpp(file = "Score_Rcpp.cpp")
# test_Score_new(x= as.matrix(p2_change),p=2,G=200, Phi=p2_a, eps= as.matrix(p2_eps) )
# #plot(msub)
#
# library(microbenchmark)
# library(ggplot2)
# mb <- microbenchmark( msub <- mosum_sub(x=p2_change,p=2,G=200, estim = "DiagC", kap = 0.3))
# mbScore <- microbenchmark( msubScore <-  mosum_sub(x=dp2_change,p=2,G=200, method="Score", estim = "DiagC", kap = 0.3))
# autoplot(mb)
# autoplot(mbScore)

##################################################


MBS_RECUR <- function(x, p, d, s, e, D, G, estim = "DiagC", var_estim = "Local", cps, stat = list(), nu=0.25, iter =1){
  if(e-s>2*G + p && iter < length(stat)){ ## segment is long enough, and recursion is shallow
    iter <- iter +1
    mod <- ar.ols(x, aic=F, order.max = p)
    mod_a <- mod$x.intercept
    eps <- mod$resid; eps[1:p,] <- 1e-4 ##solve NA
    if(p==1) mod_a <- cbind(mod_a,  matrix( mod$ar, nrow=d, ncol=d))
    if(p>1){
      for (jj in 1:p){ #collect parameters into mat
        mod_a <- cbind(mod_a,  mod$ar[jj,,])
      }
    }
    sub_vector <- get_T_RCPP(x[s:e,], p,G-p-1,mod_a,eps,estim,var_estim) ##calculate statistic on subsample ##index here
    stat[[iter]][(s+G):(e-G)] <- sub_vector[G:(e-s-G)]
    statW <- max(sub_vector) ##collect into output vector ##[which(ind==k)]
    if(statW> D){ #test
      cps_se <- get_local_maxima(sub_vector, D, G, nu=.25) #MOSUM locate change points in interval
      cps_se <- c(s, s+cps_se, e) #add start and end points
      for(jj in 2:length(cps_se)){ #BS on each interval
        cps_jj <- MBS_RECUR(x,p, d, s= cps_se[jj-1], e= cps_se[jj], D, G, estim,var_estim, cps_se, stat = stat, iter ) ##RECURSION
        stat <- cps_jj$stat
        cps_se <- union(cps_se, cps_jj$cps) ##add change points
      }
      cps <- union(cps,cps_se) ##add change points to output
    }
    #else break
  }
  return(list(cps = sort(cps), stat=stat))
}

#MBS_RECUR(dp2_change, p=2, d=3, s=100,e=400,D=2, G=100, cps = c(), stat = list(rep(0,2000),rep(0,2000),rep(0,2000)) )

MOSUMBS <- function(x, p, G, estim = "DiagC", varEstim = "Local",  alpha = 0.05){
  n <- dim(x)[1]
  d <- dim(x)[2]
  dim_warning(n,G,d,p, "Score")
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d*(d*p+1)/2 * log(log(n/G)) - log(2/3 * gamma(d*(d*p+1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  max_iter <- 3
  stat <- list(1:max_iter) #initialise statistic vector
  for (ii in 1:max_iter) {
    stat[[ii]] <- rep(0,n)
  }
  ##
  callBS <- MBS_RECUR(x,p,d,s=1,e=n,D=D_n,G,cps = c(), stat=stat, iter=0)
  cps <- callBS$cps[-c(1, length(callBS$cps))]
  #par(mfrow = c(max_iter,1))
  plot.ts(callBS$stat[[1]], ylab="Statistic") # plot test statistic
  for (ii in 2:max_iter) {
    lines(callBS$stat[[ii]], col = ii+1)
  }
  abline(h = D_n, col = "blue"); legend("topright", legend =1:max_iter,title="Iteration", fill = c(1,2:max_iter +1) ) #add threshold
  if(length(cps) > 0) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  return(list(cps = cps, stat = callBS$stat, plot=pl))
}
#msbs_test <- MOSUMBS(dp2_change, p=2, G=200)
