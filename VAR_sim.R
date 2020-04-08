## Generate VAR data with change
## parameter and error terms used throughout
a1 <- matrix(c(0.6,-0.1,-0.1,0.6),nrow=2) #regime 1
e1 <- matrix(rnorm(1000),ncol=2)
a2 <- matrix(c(-0.2,0.1,0.1,-0.2),nrow=2) #regime 2
e2 <- matrix(rnorm(1000),ncol=2)

a1 <- diag(0.7, nrow = 5, ncol = 5 ) + matrix(rnorm(5*5, 0, 0.1),ncol=5)  #regime 1 
e1 <- matrix(rnorm(5 * 1000, 0, 0.2),ncol=5)
a2 <- diag(0.2, nrow = 5, ncol = 5 ) + matrix(rnorm(5*5, 0, 0.1),ncol=5) #regime 2
e2 <- matrix(rnorm(5 * 1000, 0, 0.2),ncol=5)


## Let's start with the R version
rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(ts(simdata))
}


