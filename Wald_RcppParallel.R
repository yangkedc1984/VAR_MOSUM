library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
library(Matrix)
## RcppParallel implementation

sourceCpp(file = "Wald_RcppParallel.cpp")

ap2 <- make_a_lu(p2_change, p=2, l= 10, u= 100);ap1 <- make_a_lu(p2_change, p=1, l= 10, u= 100)

getH_ik_Wald_RCPP(as.matrix(p2_change), i=1, k=10,p=2, a = as.vector(ap2) )

makeH_k_Wald_RCPP(p2_change, k=10,p=2, a = ap2)

H_l_ex <- makeH_l_u_RCPP(p2_change, p=2, l=69, u=99, a = ap2);H_l_1 <- makeH_l_u_RCPP(p2_change, p=1, l=69, u=99, a = ap1) 
H_u_ex <- makeH_l_u_RCPP(p2_change, p=2, l=100, u=130, a = ap2);H_u_1 <- makeH_l_u_RCPP(p2_change, p=1, l=100, u=130, a = ap1) 


library(Matrix)
data <- list(
  a=matrix(1.0, 1000, 1000),
  b=matrix(2.0, 1000, 1000),
  c=matrix(3.0, 1000, 1000),
  d=matrix(4.0, 1000, 1000)
)
# z <- write_rows2(data, sapply(data, class), 1000, 1000)

get_a_lu_i_RCPP(p2_change, i=1,p=1, 10, 100) ;get_a_lu_i_RCPP(p2_change, i=1,p=2, 10, 100) ;get_a_lu_i_RCPP(p2_change, i=1,p=3, 10, 100) ;


make_a_lu_RCPP(p2_change, p=2, 10, 100)

#blockDiag(data)
V <- get_V_nk_RCPP(p2_change, p=2, 10, 100); Vbase <- get_V_nk(p2_change, p=2, 10, 100); get_V_nk(p2_change, p=4, 10, 100)

getsigma_i_kLOCAL1_RCPP(x = p2_change, i=1, k = 100, G= 30, p =2, 
                        a_upper = get_a_lu_i(p2_change, i=1, p=2, l= 100, u= 130), a_lower = get_a_lu_i(p2_change, i=1,p=2, l= 69, u= 99))
sigma_d_ex_4 <- getsigma_d_kLOCAL1_RCPP(x = p2_change, k = 100, G= 30, p =4, 
                        a_upper = make_a_lu(p2_change, p=4, l= 100, u= 130), a_lower = make_a_lu(p2_change, p=4, l= 69, u= 99))
sigma_d_ex <- getsigma_d_kLOCAL1_RCPP(x = p2_change, k = 100, G= 30, p =2, 
                        a_upper = make_a_lu(p2_change, p=2, l= 100, u= 130), a_lower = make_a_lu(p2_change, p=2, l= 69, u= 99))

#for (ii in 1:5) {
  d <- 5; p <-2
  print( (ii-1)*(d*p+1)+1-1)
  print((ii-1)*(d*p+1)+ d*p +1-1)
}

get_DiagH_Wald_RCPP(p2_change, G=30, p=2, H_l_ex, H_u_ex)

get_FullH_Wald_RCPP(p2_change, G=30, p=2, H_l_ex, H_u_ex); get_FullH_Wald_RCPP(p2_change, G=30, p=1, H_l_1, H_u_1)

get_DiagC_Wald_RCPP(x=p2_change,p=2,sigma_d = sigma_d_ex,k=100,G=30); get_DiagC_Wald_RCPP(x=p2_change,p=4,sigma_d = sigma_d_ex_4,k=100,G=30)

## evaluate Wkn (test statistic) for time k
get_Wkn_RCPP(x=p2_change,p=4,k=34,G=30, estim = "DiagC")

W_ex <- get_W_RCPP(x=p2_change,p=4,G=30, estim = "DiagC")

## TEST ----------------------------------------------------
##Wald-type test
test_Wald_new <- function(x, p, G, alpha = 0.05, estim="DiagH", ncores =1){ 
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
  Wn <- ts(get_W_RCPP(x,p,G,estim) ) #evaluate statistic at each time k
  test_stat <- max(Wn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Wn,D_n,G, nu=1/4)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  } 
  ##Plot------------------------------------
  plot(Wn, ylab="Wn") # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = Wn, cps = cps, plot = pl, estim=estim)
  return(out)
}
tW_RCPP <- test_Wald_RCPP(x=p2_change,p=2,G=200, alpha = 0.05,  estim = "DiagC") #ncores = 1,)
test_Wald_RCPP(x=p2_change,p=2,G=200, alpha = 0.05,  estim = "DiagH")
t10 <- test_Wald_RCPP(x=d10_data,p=1,G=200, alpha = 0.05,  estim = "DiagC")  #BIG





## Simulations ------
sourceCpp(file = "Wald_RcppParallel.cpp")
get_cps_RCPP(W_ex, 4.838521, 200)
test_Wald_RCPP(x=p2_change,p=2,G=200, alpha = 0.05,  estim = "DiagC")

pars <- list(a1,a2)
series_ex <- sim_data_RCPP(pars); plot.ts(series_ex)


pars <- list(a1,a2,a2,a1,a1+a2,a1%*%a2)
runsims <- var_simulate_RCPP(pars,reps=10)
