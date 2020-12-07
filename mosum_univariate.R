
mosum_univ <- function(x, p, G, method = "Wald", estim = "DiagC", varEstim = "Local",  alpha = 0.05, criterion="eps", nu=.25, rm_cross_terms =F){
  n <- dim(x)[1]
  d <- dim(x)[2] 
  
  if(rm_cross_terms) x <- remove_cross_terms(x, p, d)
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (d*p+1)/2 * log(log(n/G)) - log(2/3 * gamma( (d*p+1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector

  if(method == "Wald"){
    for (i in 1:d) {
      stat <- stat +  get_W_RCPP( as.matrix(x[,i]),p,G,estim)#W[,i] <- 
    }
  }
  if(method == "Score"){
    for(i in 1:d){
      mod <- ar(x[,i], order.max = p, demean = T, method = "ols", aic = F)
      Phi <- mod$x.intercept
      eps <- as.matrix(mod$resid); eps[1:p] <- 1e-4 ##solve NA
      if(p==1) Phi <- matrix(c(Phi,mod$ar[,,1]),1,2 ) #cbind(Phi,  matrix( mod$ar, nrow=d, ncol=d))
      if(p>1){ 
        for (jj in 1:p){ #collect parameters into mat
          Phi <- cbind(Phi,  mod$ar[jj,,])
        } 
      }
      stat <- stat +  get_T_RCPP(as.matrix(x[,i]),p,G, Phi, eps, estim)
    }  
  }
  
  cps <- c() #empty changepoint vector
  if(max(stat) > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(stat,D_n,G, nu=nu, criterion)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  } 
  ##Plot------------------------------------
  plot.ts(stat, ylab = "Statistic") # plot test statistic
  abline(h = D_n, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = D_n, mosum = stat, cps = cps, plot = pl, estim = estim)
  return(out)
}
mosum_univ(dp2_change, p=2, G=200)
mosum_univ(dp2_change, p=2, G=200, method = "Score")
mosum_univ(dp2_change, p=2, G=200, method = "Score", rm_cross_terms = T)


remove_cross_terms <- function(x,p,d){
  mod <- ar(x, order.max = p, demean = T, method = "ols", aic = F)
  
  xlist <- as.list(1:p) 
  for (jj in 1:p){ #
    aa <- matrix(mod$ar[jj,,], d,d)
    diag(aa) <- rep(0, d)
    xt <-  x %*% t(aa)
    xt <- xt[-(1:jj),]
    xt <- rbind(xt, matrix(0,jj,d))
    xlist[[jj]] <- xt
  } 
  for (jj in 1:p) {
    x <- x - xlist[[jj]]
  }
  return(x)
  
}
remove_cross_terms(dp2_change, 4, 3)
