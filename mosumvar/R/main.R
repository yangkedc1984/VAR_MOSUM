
###################
# WARNING MESSAGE #
###################

dim_warning <- function(n, G, d, p, method) {
  dim <- d*(d*p + 1 ) + d*(d+1)/2
  dimScore <- d*(d*p + 1 )/2 + d*(d+1)/2
  fl <- floor(n^(2/3))
  
  
  W3dim <- paste0("Bandwidth too small relative to model dimensions: set G > d(dp + 1) * log(d(dp + 1)) = ", d*(d*p + 1)* log(d*(d*p + 1)),
                  "\n")
  Wfl <- paste0("Bandwidth small relative to sample size: consider setting G > floor(n^(2/3)) = ", fl, "\n" )
  Wlarge <- "Large dimensions: consider `option = univariate`\n"
  
  
  if(G < dim & method == "Wald") warning(paste0("Not enough degrees of freedom for Wald method: set G > d(dp + 1) + d(d+1)/2 = ", d*(d*p + 1)* log(d*(d*p + 1)), "\n"))
  if(G < dimScore & method == "Score")warning(paste0("Not enough degrees of freedom for Score method: set G > d(dp + 1)/2 + d(d+1)/2 = ", dimScore, "\n"))
  if(G < d*(d*p + 1)* log(d*(d*p + 1)) ) warning(W3dim)
  if(G < fl ) warning(Wfl)
  if(d*(d*p + 1 ) > 30) warning(Wlarge)
  
}
# dim_warning(100, 10, 5,5)
# dim_warning(100, 145, 5,5)

#' MOSUM procedure for multiple time series
#'
#' @param x data matrix
#' @param p integer VAR model order
#' @param G integer MOSUM bandwidth
#' @param method string indicating which of `Wald` or `Score` to use
#' @param estim string estimation method
#' @param varEstim string variance estimation method
#' @param alpha Numeric significance level
#' @param criterion string location procedure
#' @param nu Numeric location procedure hyperparameter
#' @return list containing Boolean test outcome `Reject`, Numeric rejection threshold `Threshold`, 
#'  Numeric vector of test statistic `mosum`, Integer vector of estimated changepoints `cps`, Plot `plot`, 
#'  String of input estimator `estim`
#' @examples
#' data(voldata)
#' mosumvar(voldata[,2:5], 1, 250)
mosumvar <- function(x, p, G, method = c("Wald","Score")[1], estim = c("DiagC","DiagH")[1], varEstim = c("Local","Global")[1],  
                     alpha = 0.05, criterion= c("eps","eta")[1], nu=.25){
  x <- as.matrix(x)
  p <- as.integer(p)
  out <- NULL
  if(method== "Wald"){
    out <- test_Wald_new(x, p, G, alpha = alpha, estim= estim, criterion = criterion, nu=nu)
  }
  if(method== "Score"){
    out <- test_Score_new(x, p, G, alpha = alpha, estim= estim, criterion = criterion, nu=nu)
  }
  
  return(out)
}


#' Simulate multiple time series from a VAR model
#'
#' @param n integer data length
#' @param mu Numeric vector of means, defaults to zero 
#' @param Sigma error covariance matrix, defaults to identity
#' @param coeffs list or matrix of VAR coefficients; dimension `d` and lag `p` are inferred from this
#' @param error_dist string of error distribution
#' @param P1 Matrix for BEKK
#' @param Q1 Matrix for BEKK
#' @param df Integer t-distribution degrees of freedom
#' @return data frame of time series
#' @examples
#' VAR.sim(1000)
VAR.sim <- function(n, mu = NULL, Sigma = NULL, coeffs, error_dist = c("normal","t","garch")[1], P1 = NULL, Q1 = NULL, df = 1){
  if(is.matrix(coeffs)) coeffs <- list(coeffs)
  d <- ncol(coeffs[[1]])
  if(is.null(Sigma)) Sigma <- diag(1, d)
  if(is.null(mu)) mu <- rep(0, d)
  if(is.null(P1)) P1 <- matrix(1)
  if(is.null(Q1)) Q1 <- matrix(1)
  return( VAR_sim(n, mu, Sigma, coeffs, error_dist, P1, Q1, df) )
} 


#' Fit piecewise VAR model to data
#'
#' @param x data matrix
#' @param cps change point vector
#' @param p integer VAR model order (optional, uses AIC otherwise)
#' @param pen penalty scalarl; defaults to sSIC with exponent 1.01
#' @return list of model list and cost 
#' @examples
#' data(voldata)
#' run_mosum <- mosumvar(voldata[,2:5], 1, 250)
#' fit_out_model(voldata[,2:5],run_mosum$cps, p=1)
fit_out_model <- function(x,cps, p=NULL, pen = log(nrow(x))^1.01 ){
  n <- nrow(x)
  d <- ncol(x)
  starts <- c(0, cps); ends <- c(cps, n)
  
  q <- length(cps)
  RSS <- 0
  out <- as.list(1:(q+1) )
  for (ii in 1:(q+1)) {
    out[[ii]] <- ar.ols(x[(starts[ii]+1):ends[ii],] , order.max = p)
    #out[[ii]]$resid <- na.omit(out[[ii]]$resid) 
    #V <- out[[ii]]$var.pred
    RSS <- RSS + (ends[ii] - starts[ii]+1) *  norm( out[[ii]]$var.pred , type="F")^2 #   sum(diag( t(V) %*% V ))  
  }
  
  sSIC <- pen*q + (n/2) * log(RSS / n)
  return(list(model = out, sSIC = sSIC))
}