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
sigma_d_ex_4 <- getsigma_d_kLOCAL1_RCPP(x = p2_change, k = 35, G= 30, p =4, 
                        a_upper = make_a_lu(p2_change, p=4, l= 100, u= 130), a_lower = make_a_lu(p2_change, p=4, l= 69, u= 99))
sigma_d_ex <- getsigma_d_kLOCAL1_RCPP(x = p2_change, k = 100, G= 30, p =2, 
                        a_upper = make_a_lu(p2_change, p=2, l= 100, u= 130), a_lower = make_a_lu(p2_change, p=2, l= 69, u= 99))



get_DiagH_Wald_RCPP(p2_change, G=30, p=2, H_l_ex, H_u_ex)

get_FullH_Wald_RCPP(p2_change, G=30, p=2, H_l_ex, H_u_ex); get_FullH_Wald_RCPP(p2_change, G=30, p=1, H_l_1, H_u_1)

get_DiagC_Wald_RCPP(x=p2_change,p=2,sigma_d = sigma_d_ex,k=100,G=30); get_DiagC_Wald_RCPP(x=p2_change,p=4,sigma_d = sigma_d_ex_4,k=100,G=30)

## evaluate Wkn (test statistic) for time k
get_Wkn_RCPP(x=p2_change,p=4,k=100,G=30, estim = "DiagC")

W_ex <- get_W_RCPP(x=p2_change,p=4,G=30, estim = "DiagC")


# get_cps <- function(Tn, D_n, G, nu = 1/4){
#   n <- length(Tn)
#   rshift <- c(Tn[-1],0); lshift <- c(0,Tn[-n]); 
#   over <- (Tn >D_n) #indices are greater than D_n?
#   v <- which(over & lshift < D_n) #lowers
#   #v <- which(v) #append 0
#   w <- which(over & rshift < D_n) #uppers
#   #w <- which(w) #append n
#   nu_remove <- which(w-v >= nu*G) #nu(epsilon) test for distance between 
#   v_nu <- v[nu_remove]; w_nu <- c(w[nu_remove],n)
#   q <- length(v_nu) #number of CPs
#   if(q>0){
#     cps <- rep(0, q)
#     for (i in 1:q) {
#       cps[i] <- v_nu[i] + which.max(Tn[ (v_nu[i]):(w_nu[i]) ] )
#     }
#   } else {
#     cps <- NULL #if fails nu-gap
#   }  
#   return(cps)
# }

get_cps <- function(stat, D_n, G, nu = 1/4, criterion = c("eps","eta","union")[1] ){
  cps1 <- cps2 <- c() ##assign empty cps
  if(criterion =="eps" | criterion == "union"){
    sub_pairs <- get_sub_pairs(stat,D_n,G,kap=1,nu=nu) #get_sub_pairs
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
  if(criterion =="eta" | criterion == "union"){
    cps2 <- get_local_maxima(stat,D_n,G,nu=2*nu)
    q <- length(cps2)
    if(q==0) Reject <- FALSE
  }
  cps <- union(cps1,cps2) ##output union
  return(cps)
}

## TEST ----------------------------------------------------
##Wald-type test
test_Wald_new <- function(x, p, G, alpha = 0.05, estim="DiagC", criterion = "eps", nu=1/4){ 
  n <- dim(x)[1] #dimensions
  d <- dim(x)[2] 
  dim_warning(n,G,d,p,"Wald")
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d*(d*p+1)/2 * log(log(n/G)) - log(2/3 * gamma(d*(d*p+1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  Wn <- get_W_RCPP(x,p,G,estim)  #evaluate statistic at each time k
  test_stat <- max(Wn)
  cps <- c() #empty changepoint vector
  if(test_stat > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Wn,D_n,G, nu=nu, criterion)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  } 
  ##Plot------------------------------------
  plot.ts(Wn, ylab="Wn") # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = Wn, cps = cps, plot = pl, estim=estim)
  return(out)
}
tW_RCPP <- test_Wald_new(x=p2_change,p=2,G=200, alpha = 0.05,  estim = "DiagC") #ncores = 1,)
test_Wald_RCPP(x=p2_change,p=2,G=200, alpha = 0.05,  estim = "DiagH")
t10 <- test_Wald_RCPP(x=d10_data,p=1,G=200, alpha = 0.05,  estim = "DiagC")  #BIG




## Simulations ------
sourceCpp(file = "Wald_RcppParallel.cpp")
get_cps_RCPP(W_ex, 4.838521, 200)
test_Wald_new(x=p2_change,p=2,G=200, alpha = 0.05,  estim = "DiagC")

pars <- list(a1,a2)
series_ex <- sim_data_RCPP(pars); plot.ts(series_ex)


pars <- list(a1,a2,a2,a1,a1+a2,a1%*%a2)
runsims <- var_simulate_RCPP(pars,reps=10)
