

#Evaluate CUSUM at b in window [s,e-1]
CUSUM <- function(x, b, s, e) {
  n <- e - s +1 #range
  Xb <- sqrt( (e-b) / (n*(b-s+1)) ) * sum(x[s:b]) - sqrt( (b-s+1) / (n*(e-b)) ) * sum(x[(b+1):e])  #calculate CUSUM
  Xb #return CUSUM
}


#Call CUSUM for each index in winddow [s, e-1]
call.CUSUM <- function(x, s, e){
  B <- matrix(s:(e-1)) #index for window [s,e-1]
  XB <- apply(X = B, MARGIN = 1, FUN = CUSUM, x=x, s=s, e=e) #find CUSUM across index
  XB
}

BINSEG <- function(x, s, e, z, cps = c()) {
  
  if(e-s<1){break}
  
  x_CUSUM <- abs(call.CUSUM(x, s, e))
  b0 <- which.max(x_CUSUM)
  
  if(x_CUSUM[b0] > z){
    append(cps, b0) #add to change points
    BINSEG(x, s, b0, z, cps)
    BINSEG(x, b0 + 1, z, cps)
  } else {break}
  
  return(cps)
}


x <- c(rep(0, 100), rep(10, 100)) + rnorm(200, 0, 10)


x_CUSUM <- call.CUSUM(x, 2, length(x) -1 )
plot.ts(x_CUSUM)

BINSEG(x, 2, length(x) - 1, 0.9)
