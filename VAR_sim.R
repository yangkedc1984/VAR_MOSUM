## Generate VAR data with change
## parameter and error terms used throughout
# a1 <- matrix(c(0.6,-0.1,-0.1,0.6),nrow=2) #regime 1
# e1 <- matrix(rnorm(1000),ncol=2)
# a2 <- matrix(c(-0.2,0.1,0.1,-0.2),nrow=2) #regime 2
# e2 <- matrix(rnorm(1000),ncol=2)
# 
# ##change
# a1 <- diag(0.6, nrow = 5, ncol = 5 ) + 0.05  #regime 1 
# e1 <- matrix(rnorm(5 * 1000, 0, 0.2),ncol=5)
# a2 <- diag(-0.4, nrow = 5, ncol = 5 ) - 0.03 #regime 2
# e2 <- matrix(rnorm(5 * 1000, 0, 0.2),ncol=5)
# 
# 
# #univ 
# a1_univ <- 0.7; a2_univ <- -0.7
# e1_univ <- rnorm(500,0,0.07); e2_univ <- rnorm(500,0,0.07)

## Let's start with the R version
rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  simdata[1,] <- errors[1,]
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(ts(simdata))
}
## p=2
rSim_p2 <- function(coeff, coeff2,  errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  simdata[1:2,] <- errors[1:2,]
  for (row in 3:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + coeff2 %*% simdata[(row-2),] + errors[row,]
  }
  return(ts(simdata))
}
#rSim_p2(a1,a2,e1)

## p=5
rSim_p <- function(coeff_list,  errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  p <- length(coeff_list)
  simdata[1:p,] <- errors[1:p,]
  for (row in (p+1):nrow(errors)) {
    simdata[row,] <-  errors[row,]
    for (i in 1:p) {
      simdata[row,] <-simdata[row,] + coeff_list[[i]] %*% simdata[(row-i),] 
    }
  }
  return(ts(simdata))
}


Alist <- list( matrix( c(0.7,-.1,.7,-.1), nrow = 2, ncol = 2),  matrix( -c(0.1,.1,.1,.1), nrow = 2, ncol = 2) )
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
library(Matrix)
sourceCpp(file = "VAR_sim.cpp")
##

plot.ts(
  VAR_sim(n=1000, mu = c(0,0), Sigma = diag(1,2,2), coeffs = list( diag(0.9,2,2) ), error_dist = "normal", P1= diag(1), Q1 = diag(1) ) )

plot.ts(
  VAR_sim(n=1000, mu = c(0,0), Sigma = diag(0.5,2,2), coeffs = list( diag(0.9,2,2) ), error_dist = "t", df = 100 ) )

plot.ts(
  VAR_sim(n=1000, mu = c(0,0), Sigma = diag(0.5,2,2), coeffs = list( diag(0.9,2,2) ), error_dist = "garch", df = 100,
          P1= diag(0.5, 2,2), Q1 = diag(0.5, 2,2) ) )

### mb


#library(microbenchmark)
#microbenchmark()