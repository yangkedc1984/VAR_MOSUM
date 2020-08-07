library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
library(Matrix)
## RcppParallel implementation

## Sim ------------------
rData1 <- rSim(a1, e1)
rData1[1,] <- runif(5, 0, 0.02) #prevent NaN
rData2 <- rSim(a2, e2)
rData2[1,] <- runif(5, 0, 0.02)

var_change <- ts(rbind(rData1+ 0, rData2- 0) )
change_model <- ar(var_change, order.max = 1, demean = T, method = "ols")
a_change <- cbind(change_model$x.intercept, matrix( change_model$ar, nrow=5, ncol=5))
eps_change <- change_model$resid
plot(var_change)
##change p=2 data------------------------- 
#
p2_Data1 <- rSim_p2(a1, a2, e1)
p2_Data2 <- rSim_p2(a2, -a1, e2)

p2_change <- ts(rbind(p2_Data1, p2_Data2) )
p2_model <- ar(p2_change, order.max = 2, demean = T, method = "ols")
p2_a <- cbind(p2_model$x.intercept,  p2_model$ar[1,,],  p2_model$ar[2,,])
p2_eps <- p2_model$resid
plot(p2_change)
#------------------------------------
sourceCpp(file = "Score_Rcpp.cpp")
getH_ik(var_change,i=1,k=100,p=1, Phi=a_change, eps=eps_change)

makeH_k_Wald_RCPP(var_change, k=100,p=1, Phi=a_change, eps=eps_change)

hall <- makeH_all_RCPP(var_change, p=1, G=200, Phi=a_change, eps=eps_change)

get_DiagH_RCPP(var_change, k=100,G=30, p=1, hall)

get_FullH_RCPP(var_change, k=100,G=30, hall)

getA_RCPP(var_change, k=300, G=200, p=1, eps_change, hall) ; getA(var_change, k=300, G=200, p=1, a_change, eps_change, hall)

getsigma_iGlobal_RCPP(eps_change,1,1)
getsigma_dGlobal_RCPP(eps_change,1)
getsigma_iLocal_RCPP(eps_change,1,31,30)
sdl <- getsigma_dLocal_RCPP(eps_change, 1, 200)

get_DiagC_RCPP(var_change, p=1, sdl[300,], k=300, G=200 )

tkn <- get_Tkn_RCPP(var_change,  k=300,p=1,G=200,Phi=a_change, eps=eps_change,h_all=hall, as.matrix(sdl[300,]), estim = "DiagC", var_estim = "Global"  )
tn <- get_T_RCPP(var_change,p=1,G=200,Phi=a_change, eps=eps_change)
tn <- get_T_RCPP(var_change,p=1,G=200,Phi=a_change, eps=eps_change, estim = "DiagC", var_estim = "Global")
tn <- get_T_RCPP(p2_change,p=2,G=200,Phi=p2_a, eps=p2_eps, estim = "DiagH", var_estim = "Global")
ts.plot(tn)

##Score-type test
test_Score_new <- function(x, p, G, Phi, eps, alpha = 0.05, estim="DiagC",var_estim = "Local"){ 
  if(is.null(dim(x)) || dim(x)[2] == 1 ) {x <- matrix(x); Phi <- matrix(Phi); eps <- matrix(eps)} #handle univariate case
  n <- dim(x)[1] #dimensions
  d <- dim(x)[2] 
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d*(d*p+1)/2 * log(log(n/G)) - log(2/3 * gamma(d*(d*p+1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  Tn <- ts(get_T_RCPP(x,p,G,Phi,eps,estim,var_estim)) #evaluate statistic at each time k
  test_stat <- max(Tn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Tn,D_n,G, nu=1/4)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
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
