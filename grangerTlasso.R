#################################################################################
#            SOME MATRIX MANIPULATION FUNCTIONS
#################################################################################
adj_right2left = function(adjmat){
  ntime = dim(adjmat)[1]
  p = dim(adjmat)[2]
  adjout = array(0, c(p, p, ntime))
  for (tt in 1:ntime)
    for (i in 1:p)
      for (j in 1:p)
        adjout[i,j,ntime-tt+1] = adjmat[tt, i,j]  #p-j+1
  
  return(adjout)
}


adj_left2right = function(adjmat){
  ntime = dim(adjmat)[3]
  p = dim(adjmat)[2]
  adjout = array(0, c(ntime, p, p))
  for (tt in 1:ntime)
    for (i in 1:p)
      for (j in 1:p)
        adjout[tt, i,j] = adjmat[i,j,ntime-tt+1]   #p-j+1
  
  return(adjout)
}

MatImage <- function(AA, maintitle=NULL, margin=c(2,1,2,1)) {
  d1 = dim(AA)[1]
  d2 = dim(AA)[2]
  
  par(mar=margin)
  image( t(AA[d1:1,]) , col=gray(100:0/100), axes=FALSE )
  box()
  title(maintitle, cex.main = 2,   font.main= 1)
  
  axis(2, at=c(d1,1)/d1, label=c(d1,1))
  
  if (d1 < d2){
    axis(1, at=seq(d2,1,by=-d1)/d2, label=seq(1,d2,by=d1))
  }else{
    axis(1, at=c(d2,1)/d2, label=c(1,d2))
  }
}

array2mat <-
  function(
    myArray,		#input array
    xtnd=2		#dimension to extend the array over
  ){
    ####### START OF FN #######
    d1 = dim(myArray)[1]
    d2 = dim(myArray)[2]
    d3 = dim(myArray)[3]
    
    if (is.na(d3)){
      myMat = myArray
      return(myMat)
    }
    
    if (xtnd==1){
      myMat = matrix(0,d1*d3,d2)
    }else{
      myMat = matrix(0,d1,d2*d3)
    }
    
    for(i in 1:d3){
      low1 = (i-1)*d1+1
      high1 = i*d1
      low2 = (i-1)*d2+1
      high2 = i*d2
      
      if (xtnd==1){
        myMat[low1:high1,] = myArray[,,i]
      }else{
        myMat[,low2:high2] = myArray[,,i]
      }
    }
    
    return(myMat)
  }


#' Use truncating lasso to estimate graphical Granger causality https://rdrr.io/github/njetzel/ngc/src/R/grangerTlasso.R
#' @param X input array
#' @param d number of time lags to consider
#' @param group group indices
#' @param typeIerr acceptable type I error rate
#' @param typeIIerr acceptable type II error rate
#' @param weights weights for adaptive lasso
#' @param Alasso.power power for adaptive lasso
#' @param eps epsilon used for truncating penalty
#' @param tol tolerance used for BR algorithm
#' @param thresh values smaller than thresh are set to zero
#' @param initialEst initial estimate for BR alg
#' @param maxIter maximum number of iterations for BR alg
#' @return a list including a matrix of estimated coefficients, final lambda value, and time series order estimate
grangerTlasso <-
  function(
    X, 					#input array dim=(n,p,T), last time=Y
    d = NULL, 		#number of time lags to consider
    group = NULL,
    typeIerr = 0.10, 			#sig level for lambda (...Method=errBased)
    typeIIerr = 0.10,			#acceptable type II error
    weights = NULL,  #matrix of weights for Alasso
    Alasso.power = 1,			#power for Adaptive lasso
    eps = 1e-8,				#epsilon, used for truncating penalty
    tol = 1e-2,				#tolerance, for BR alg
    thresh = 1e-4,		#values smaller than thresh are set to zero
    initialEst=NULL, 	#initial estimate for BR alg
    maxIter=100				#maximum iteration for the BR alg
  ){
    ####### START OF FN #######
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]
    
    useAlasso <- !is.null(weights)
    
    Onep <- matrix(1,p,p)
    
    newEst <- initialEst
    if (is.null(newEst)){
      newEst <-  array( 0, c(p, p, d) )
    }
    
    ##scale the X matrix
    for (i in 1:(tp-1)){
      X[,,i] <- scale( X[,,i] )#*sqrt(n/(n-1))
    }
    
    Y <- X[,,tp]		#for compatiblity with previous codes
    YY <- Y
    
    diff <- tol + 1
    iter <- 0			#no of iterations in the while loop
    FLAG <- TRUE		#flag for whether to stop the while loop
    d <- d		#no of lags to consider in estimation
    lambda <- NULL
    #while loop for BR algorithm
    #cat("BR alg started! \n")
    while ( (diff > tol) && (iter<=maxIter) && (FLAG==TRUE) ){
      
      oldEst <- newEst
      iter <- iter + 1
      #cat("BR loop No: ", iter, "\n")
      
      ##estimation loop
      jj <- 1
      CONTINUE <- TRUE		#flag for continuation of inner estimation loop
      while (jj <= d){
        XX <- X[,,(tp-jj)]
        ##Find the residual as Y (similar to Gauss-Sidel)
        theta <- array2mat(newEst[,,-(tp-jj)])
        for (i in 1:p){
          YY[,i] <- Y[,i] -
            array2mat(X[, , -c(tp-jj,tp)]) %*% theta[i,]
        }
        
        ## calculate truncating factor
        if (jj>1){
          previousNNZ <- sum( abs(newEst[,,(tp-jj+1)]) > thresh )/(p^2)
          if( previousNNZ < (typeIIerr/(tp-jj)) ){
            
            newEst[,,1:(tp-jj)] <- matrix(0,p,p*(tp-jj))
            d <- jj - 1
            CONTINUE <- FALSE
          }
        }
        
        ## calculate est
        if (CONTINUE==TRUE){
          if (!useAlasso){
            tempp <- pldag.set(X1=XX, X2=YY, group = group,
                               sigLevel=typeIerr, useWghts=FALSE,
                               wantScale=TRUE, useTLASSO=TRUE, d=d)
            
            newEst[,,(tp-jj)] <- tempp$AA
          }else{
            
            W <- (abs(weights[,,(tp-jj)]) + eps*Onep)^(-Alasso.power)
            if (sum(W < 1) > 0){
              W[(W < 1)] <- 1
            }
            tempp <-pldag.set(XX, YY, group = group,
                              sigLevel=typeIerr, useWghts=TRUE,
                              wghts=W, wantScale=TRUE, useTLASSO=TRUE, d=d)
            
            newEst[,,(tp-jj)] <- tempp$AA
          }#endif(!useAlasso)
        }#endif(CONTINUE)
        jj <- jj + 1
        
      }#endwhile(jj)
      
      diff <- sum( abs(newEst-oldEst) ) / sum( abs(newEst)>eps )
      if ((diff == Inf) || is.na(diff))
      {
        diff = tol + 1
      }
    }#endwhile()
    lambda <- tempp$lambda
    intercepts <- tempp$intercepts
    
    if (diff > tol){
      cat("WARNING: BR Alg not converged! Increase maxIter. \n")
    }
    
    return(list(estMat=newEst, lambda=lambda, tsOrder=jj-1, intercepts = intercepts))
  }


grangerTlasso(graphtensor, d=5)
