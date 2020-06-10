## data generation

x <- rbind(cbind(arima.sim(list(ar = .5), n = 1000), arima.sim(list(ar = -.1), n = 1000), arima.sim(list(ar = .5), n = 1000)),
           cbind(arima.sim(list(ar = -.5), n = 1000), arima.sim(list(ar = -.7), n = 1000), arima.sim(list(ar = .1), n = 1000)),
           cbind(arima.sim(list(ar = .1), n = 1000), arima.sim(list(ar = .1), n = 1000), arima.sim(list(ar = .5), n = 1000)))

x <- rbind(cbind(arima.sim(list(ar = .5), n = 3000), arima.sim(list(ar = -.1), n = 3000), arima.sim(list(ar = .5), n = 3000)))

# 5.1.2

n <- 2000
d <- 3
burnin <- 100
cp <- 500 * (1:3)
brks <- c(0, cp + burnin, n + burnin)
L <- c(list(c(list(matrix(c(.5, .1, .1, .1, .5, .1, .1, .1, .5), byrow = TRUE, nrow = 3)), list(matrix(c(-1, 1, 1, 1, -1, 1, 1, 1, -1) * .2, byrow = TRUE, nrow = 3)))),
       list(c(list(matrix(c(-1, 1, 1, 1, -1, 1, 1, 1, -1) * .2, byrow = TRUE, nrow = 3)), list(matrix(c(.5, .1, .1, .1, .5, .1, .1, .1, .5), byrow = TRUE, nrow = 3)))),
       list(c(list(matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1) * .4, byrow = TRUE, nrow = 3)), list(diag(rep(-.2, 3))))),
       list(c(list(matrix(c(.5, .1, .1, .1, .5, .1, .1, .1, .5), byrow = TRUE, nrow = 3)), list(matrix(c(-1, 1, 1, 1, -1, 1, 1, 1, -1) * .2, byrow = TRUE, nrow = 3)))))

x <- matrix(rnorm((n + burnin) * d), ncol = d)
for(jj in 2:(length(cp) + 1)){
  int <- setdiff((brks[jj] + 1):brks[jj + 1], 1:2)
  for(tt in int) x[tt, ] <- L[[jj]][[1]] %*% x[tt - 1, ] + L[[jj]][[2]] %*% x[tt - 2, ] + x[tt, ]
}
x <- x[-(1:burnin), ]


# 5.1.3

n <- 2000
d <- 2
burnin <- 100
cp <- 500 * (1:3)
brks <- c(0, cp + burnin, n + burnin)
L <- c(list(matrix(c(-1, -1, 1, 1) * .75, byrow = TRUE, nrow = 2)), 
       list(matrix(-c(-1, -1, 1, 1) * .25, byrow = TRUE, nrow = 2)),
       list(matrix(c(-1, -1, 1, 1) * .25, byrow = TRUE, nrow = 2)),
       list(matrix(-c(-1, -1, 1, 1) * .75, byrow = TRUE, nrow = 2)))

x <- matrix(rnorm((n + burnin) * d), ncol = d)
for(jj in 1:(length(cp) + 1)){
  int <- setdiff((brks[jj] + 1):brks[jj + 1], 1)
  for(tt in int) x[tt, ] <- L[[jj]] %*% x[tt - 1, ] + x[tt, ]
}
x <- x[-(1:burnin), ]

##### some explanatory analysis, plotting, ar fitting, etc.

ts.plot(x[, 2])
ts.plot(x[, 1] * x[, 2])

colnames(x) <- c('ts1', 'ts2', 'ts3')[1:d]

# matplot(x, type = 'l')

fit <- ar(x, demean = TRUE, method = 'ols')
fit$order




p <- 1
fit <- ar(x,aic = F,  order = p, demean = TRUE, method = 'ols')

# z <- fit$resid[-(1:p), ]
fit$ar

A <-  cbind(fit$x.intercept, matrix(fit$ar,2,2))
  

###### start

p <- 1
n <- dim(x)[1]; d <- dim(x)[2]

G <- 300

alpha <- max(.01,  (1 - (1/exp(exp(-sqrt(2*log(n/G)) * sqrt(G))))^2))

D <- d*(d*p + 1)
c_alpha <- -log(log(1/sqrt(1 - alpha))) 
a <- sqrt(2*log(n/G)) #test transform multipliers
b <- 2*log(n/G) + D/2 * log(log(n/G)) - log(2/3 * gamma(D/2))

print(crit.val <- max(sqrt(2 * log(n)), (b + c_alpha)/a))

X <- c()
for(ii in (p + 1):n) X <- rbind(X, c(1, x[ii - (1:p), ]))
y <- x[-(1:p), ]

solve(t(X)%*%X)%*%t(X)%*%y

z <- y - X %*% solve(t(X)%*%X)%*%t(X)%*%y

H <- array(0, dim = c(n - p, dim(X)[2], d))
for(jj in 1:d) H[,, jj] <- X * z[, jj]

#### Sigma (obtained locally) * C (globally/locally, determined by do.local)

do.local <- !TRUE

tmp <- stat <- rep(0, n - p)
ii <- G
int <- (ii - G + 1):(ii + G)
left <- apply(H[1:G,, ], c(2, 3), sum)
right <- apply(H[G + 1:G,, ], c(2, 3), sum)
A <- (left - right)/sqrt(2 * G)
sigma <- diag(t(z[int, ]) %*% z[int, ]/(2 * G))
if(do.local) C <- t(X[int, ]) %*% X[int, ] / (2 * G) else C <- t(X) %*% X / (n - p) 
Cinv <- solve(C)
for(jj in 1:d) stat[ii] <- stat[ii] + sum(A[, jj] * t(A[, jj] * Cinv))/sigma[jj]

tmp[ii] <- sum(A^2)

for(ii in (G + 1):(n - p - G)){
  int <- (ii - G + 1):(ii + G)
  left <- left - H[ii - G,,] + H[ii,,]
  right <- right - H[ii,,] + H[ii + G,,]
  A <- (left - right)/sqrt(2 * G)
  sigma <- diag(cov(z[int, ]))
  if(do.local){
    C <- t(X[int, ]) %*% X[int, ] / (2*G)
    Cinv <- solve(C)  # can adopt the woodbury formula for efficiency here
  }
  for(jj in 1:d) stat[ii] <- stat[ii] + sum(A[, jj] * t(A[, jj] * Cinv))/sigma[jj]
  
  tmp[ii] <- sum(A^2)
}

ts.plot(sqrt(tmp))

ts.plot(sqrt(stat)); abline(v = cp, h = crit.val, col = 2)

#### H (diagonal / not diagonal, determined by do.diag)

do.diag <- !TRUE

stat <- rep(0, n - p)

int <- 1:(2 * G)
left <- apply(H[1:G,, ], c(2, 3), sum)
right <- apply(H[G + 1:G,, ], c(2, 3), sum)
A <- (left - right)/sqrt(2 * G)
# XX <- solve(t(X[1:(2*G), ]) %*% X[1:(2*G), ]) 
# for(jj in 1:d) stat[G] <- stat[G] + sum(A[, jj] * t(A[, jj] * XX))/sigma[jj]
for(jj in 1:d){
  HH <- cov(H[int[1:G],, jj]) + cov(H[int[G + 1:G],, jj])
  if(do.diag){
    HH <- diag(1/diag(HH))
  } else HH <- solve(HH) 
  stat[ii] <- stat[ii] +  sum(A[, jj] * t(A[, jj] * HH))
}

for(ii in (G + 1):(n - p - G)){
  int <- (ii - G + 1):(ii + G)
  left <- left - H[ii - G,,] + H[ii,,]
  right <- right - H[ii,,] + H[ii + G,,]
  A <- (left - right)/sqrt(2 * G)
  # XX <- solve(t(X[int, ]) %*% X[int, ]/(2 * G))
  # for(jj in 1:d) stat[ii] <- stat[ii] +  sum(A[, jj] * t(A[, jj] * XX))/sigma[jj]
  for(jj in 1:d){
    HH <- (cov(H[int[1:G],, jj]) + cov(H[int[G + 1:G],, jj]))/2
    if(do.diag){
      HH <- diag(1/diag(HH))
    } else HH <- solve(HH) 
    stat[ii] <- stat[ii] +  sum(A[, jj] * t(A[, jj] * HH))
  }
}

ts.plot(sqrt(stat)); abline(v = cp, h = crit.val, col = 2)

## IGNORE
# 
# jj <- 1
# 
# ii <- 1000
# int <- (ii - G + 1):(ii + G)
# sigma <- diag(t(z[int, ]) %*% z[int, ])/(2 * G)
# if(do.local){
#   C <- t(X[int, ]) %*% X[int, ] / (2*G)
#   Cinv <- solve(C)  # can adopt the woodbury formula for efficiency here
# } else{
#   C <- t(X) %*% X / (n - p) 
#   Cinv <- solve(C)
# }
# Cinv/sigma[jj]
# 
# HH <- (cov(H[int[1:G],, jj]) + cov(H[int[G + 1:G],, jj])) * .5
# if(do.diag){
#   HH <- diag(1/diag(HH))
# } else HH <- solve(HH) 
# 
# HH

