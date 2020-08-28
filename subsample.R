## Subsampling algorithm
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
library(Matrix)
sourceCpp(file = "Wald_RcppParallel.cpp")
sourceCpp(file = "Score_Rcpp.cpp")
##

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
  over <- (Tn >D_n) #indices are greater than D_n?
  v <- which(over & lshift < D_n )#& lshift >0 ) #lowers
  #v <- which(v) #append 0
  w <- which(over & rshift < D_n )#& rshift >0) #uppers
  #w <- which(w) #append n
  nu_remove <- which(w-v >= nu*R*G) #nu(epsilon) test for distance between 
  v_nu <- v[nu_remove]; w_nu <-w[nu_remove] #c(w[nu_remove],n)
  # q <- length(v_nu) #number of CPs
  # if(q>0){
  #   cps <- rep(0, q)
  #   for (i in 1:q) {
  #     cps[i] <- v_nu[i] + which.max(Tn[ (v_nu[i]):(w_nu[i]) ] )
  #   }
  # } else {
  #   cps <- NULL #if fails nu-gap
  # }  
  sub_pairs <- cbind(v_nu,w_nu)  
  return(sub_pairs)
  #plot(lshift); lines(Tn)
}
get_sub_pairs(Tn = msub$mosum, D_n = 5, G= 200, kap = 0.3, nu = .25)

##get subsample change point estimates

mosum_sub <- function(x, p, G, estim = "DiagC", varEstim = "Local", kap = 0.1,  alpha = 0.05){
  n <- dim(x)[1]
  d <- dim(x)[2] 
  #nu <- 0.25
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
  statW <- vapply(ind, get_Wkn_RCPP, FUN.VALUE = double(1), x=x,p=p,G=G, estim=estim)#get_W_RCPP(x[ind,],p,G,estim)
  stat[ind] <- statW
  #return(stat)
  test_stat <- max(stat)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    #tlist <- list(which(test_stat> D_n))
    #for (i in 0:G ){ #collect times around significant Wk
    #  tlist <- append(tlist, c(which(test_stat> D_n) + i,which(test_stat> D_n) -i)  )
    #}
    times <- which(stat> D_n)#times <- sort(Reduce(union,tlist))
    tlist <- split(times, cumsum(c(1, diff(times) != R))) #split into list of consecutive regions
    
    for (i in 1:length(tlist)) {
    interval <- (min(tlist[[i]]) - 2*G):(max(tlist[[i]]) + 2*G)
    #fit var model
    mod <- ar(x[interval,], order.max = p, demean = T, method = "ols")
    mod_a <- mod$x.intercept
    eps <- mod$resid
     for (jj in 1:p){ #collect parameters into mat
       mod_a <- cbind(mod_a,  mod$ar[jj,,])
     } 
    #interval.sort <- sort( union(interval-p,interval))
    #stat[(min(tlist[[i]])):(max(tlist[[i]]))] <- get_T_RCPP( as.matrix(x[interval,]),p,G,Phi= as.matrix(mod_a), eps=as.matrix(eps),estim = estim) #overwrite statistic
    stat[interval[(1*G):(length(interval)-1*G )]] <- get_T_RCPP( 
      as.matrix(x[interval,]),p,G,Phi= as.matrix(mod_a), eps=as.matrix(eps),estim = estim)[(1*G):(length(interval)-1*G )] #overwrite statistic
    }
    
    sub_pairs <- get_sub_pairs(stat,D_n,G,kap=kap,nu=0.25) #get_sub_pairs
    q <- dim(sub_pairs)[1]
    if(q==0) Reject <- FALSE
    else if (q>0){ ## locate cps
     for (ii in 1:q) {
       interval <- sub_pairs[ii,1]:sub_pairs[ii,2]
       kk <- which.max(stat[interval]) #internal cp location
       cps[ii] <- kk + sub_pairs[ii,1] #- G-p
       #stat[sub_pairs[ii,1]:sub_pairs[ii,2]] <- statT[(G+p):(sub_pairs[ii,2]-sub_pairs[ii,1] +G+p )]
     }
    }
    #cps <- get_cps(Wn,D_n,G, nu=1/4)
    #if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  } 
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

msub <- mosum_sub(x=p2_change,p=2,G=200, estim = "DiagC", kap = 0.6)
mmm <- mosum_sub(x=dp2_change,p=2,G=200, estim = "DiagC", kap = 1)
  test_Wald_new(x=dp2_change,p=2,G=200, alpha = 0.05,  estim = "DiagC")
  
sourceCpp(file = "Score_Rcpp.cpp")
test_Score_new(x= as.matrix(p2_change),p=2,G=200, Phi=p2_a, eps= as.matrix(p2_eps) )
#plot(msub)

library(microbenchmark)
library(ggplot2)
mb <- microbenchmark( msub <- mosum_sub(x=p2_change,p=2,G=200, estim = "DiagC", kap = 0.3))
autoplot(mb)
