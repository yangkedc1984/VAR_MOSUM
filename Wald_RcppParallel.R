library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
## RcppParallel implementation


##estimating fn for channel i, time k
# getH_ik_Wald_RCPP <- function(x, i,k,p, a) { 
#   d <- dim(x)[2]
#   ai <- a[((i-1)*(d*p+1)+1):((i-1)*(d*p+1)+ d*p +1)] #select params for channel i #accounts for intercept
#   if(p==1)X <- c(1, x[(k-1),]) #with intercept X <- x[(k-1),]
#   if(p==2)X <- c(1, x[(k-1),], x[(k-2),])
#   if(p==3)X <- c(1, x[(k-1),], x[(k-2),],x[(k-3),])
#   if(p==4)X <- c(1, x[(k-1),], x[(k-2),],x[(k-3),],x[(k-4),])
#   y <- as.numeric(x[k,i])
#   e <- as.numeric(x[k,i] - ai %*% X) #eps[k,i] #p=1 residual
#   H_ik <- t(- y*X + ai%*%X%*%t(X) - e*X) #*2
#   return( H_ik)
# }
sourceCpp(file = "Wald_RcppParallel.cpp")

ap2 <- make_a_lu(p2_change, p=2, l= 10, u= 100);ap1 <- make_a_lu(p2_change, p=1, l= 10, u= 100)

getH_ik_Wald_RCPP(as.matrix(p2_change), i=1, k=10,p=2, a = as.vector(ap2) )

makeH_k_Wald_RCPP(p2_change, k=10,p=2, a = ap2)

H_l_u <- makeH_l_u_RCPP(p2_change, p=2, l=10, u=100, a = ap2)
H_l_u

library(Matrix)
data <- list(
  a=matrix(1.0, 1000, 1000),
  b=matrix(2.0, 1000, 1000),
  c=matrix(3.0, 1000, 1000),
  d=matrix(4.0, 1000, 1000)
)
# z <- write_rows2(data, sapply(data, class), 1000, 1000)

get_a_lu_i_RCPP(p2_change, i=1,p=1, 10, 100) ;get_a_lu_i_RCPP(p2_change, i=1,p=2, 10, 100) ;

make_a_lu_RCPP(p2_change, p=1, 10, 100)

#blockDiag(data)
V <- get_V_nk_RCPP(p2_change, p=2, 10, 100)

getsigma_i_kLOCAL1_RCPP(x = p2_change, i=1, k = 100, G= 30, p =2, 
                        a_upper = get_a_lu_i(p2_change, i=1, p=2, l= 100, u= 130), a_lower = get_a_lu_i(p2_change, i=1,p=2, l= 69, u= 99))
getsigma_d_kLOCAL1_RCPP(x = p2_change, k = 100, G= 30, p =2, 
                        a_upper = make_a_lu(p2_change, p=2, l= 100, u= 130), a_lower = make_a_lu(p2_change, p=2, l= 69, u= 99))
