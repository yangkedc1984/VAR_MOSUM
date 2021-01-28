# ## Multi scale
# library(Rcpp)
# library(RcppParallel)
# library(RcppArmadillo)
# library(Matrix)
# sourceCpp(file = "Wald_RcppParallel.cpp")
# sourceCpp(file = "Score_Rcpp.cpp")
# ###
#
# Gset <- c(50, 100, 150)

# MFA <- function(x, p, Gset, estim = "DiagC",  alpha = 0.05){
#   Gset <- sort(Gset) #into ascending order
#   Glen <- length(Gset)
#   cps <- c()
#   Reject <- F
#   tests <- list(1:Glen)
#   for (i in 1:Glen){
#     tests[[i]] <-  test_Wald_new(x, p, Gset[i], estim = estim,  alpha = alpha)
#     if(tests[[i]]$Reject) Reject <- TRUE
#   }
#   if(Reject){
#     cps <- tests[[1]]$cps
#     if(Glen > 1){
#       for (i in 2:Glen){
#         K <- tests[[i]]$cps
#         for (j in 1:length(K)) {
#           if(min(abs(K[j] - cps)) > G) append(cps, K[j])
#         }
#       }
#     }
#   }
#   plot.ts(x, ylab="Series") # plot series
#   #abline(h = D_n, col = "blue") #add threshold
#   if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
#   pl <- recordPlot()
#   return(list(Reject, cps, q= length(cps), pl=pl))
# }
# MFA(x=dp2_change, p=2,Gset=c(90,200,300), estim = "DiagC")


#MFA_RCPP(x=dp2_change, p=2,Gset=c(200,300), estim = "DiagC")

#' MOSUM procedure for multiple time series at multiple scales
#'
#' @param x data matrix
#' @param p integer VAR model order
#' @param Gset integer vector of MOSUM bandwidths
#' @param method string indicating which of `Wald` or `Score` to use
#' @param estim string estimation method
#' @param varEstim string variance estimation method
#' @param alpha Numeric significance level
#' @return list containing Boolean test outcome `Reject`, Numeric rejection threshold `Threshold`, 
#'  Numeric vector of test statistic `mosum`, Integer vector of estimated changepoints `cps`, Plot `plot`, 
#'  String of input estimator `estim`
#' @examples
#' data(voldata)
#' MFA(voldata[,2:5], 1, c(100, 250, 400) )
MFA <- function(x, p, Gset, test = c("Wald","Score")[1],  estim = c("DiagC","DiagH")[1],  alpha = 0.05){
   x <- as.matrix(x)
   p <- integer(p)
   out <- NULL
 if(test=="Wald") out <- MFA_Wald(x,p,Gset,estim,alpha)
 if(test=="Score"){
   mod <- ar(x, order.max = p, demean = T, method = "ols", aic = F)
   Phi <- mod$x.intercept
   eps <- mod$resid; eps[1:p,] <- 1e-4 ##solve NA
   if(p==1) Phi <- cbind(Phi,  matrix( mod$ar, nrow=ncol(x), ncol=ncol(x)))
   if(p>1){
     for (jj in 1:p){ #collect parameters into mat
       Phi <- cbind(Phi,  mod$ar[jj,,])
     }
   }
   out <- MFA_Score(x,p,Gset, Phi, eps, estim,alpha)
 }
 plot.ts(x)
 if(out$Reject) abline(v = out$ChangePoints, col = "red")  #if rejecting H0, add estimated cps
 out$Changepoints <- sort(out$Changepoints)
 #pl <- recordPlot()
 return(out)
}
#MFA(x=dp2_change, p=2,Gset=c(50,200,300), test = "Wald", estim = "DiagC")
