
mosum_univ <- function(x, p, G, method = "Wald", estim = "DiagC", varEstim = "Local",  alpha = 0.05, criterion="eps", nu=.25,
                       rm_cross_terms =F, do_bootstrap = F, M = 1000, global_resids = F){
  n <- dim(x)[1]
  d <- dim(x)[2] 
  if(global_resids) {
    eps_global <-  (ar(x, order.max = p, demean = T, method = "ols", aic = F)$resid)
    #varEstim <- "Global"
  }  
  if(rm_cross_terms) x <- remove_cross_terms(x, p, d)
  ##Test setup----------------------------
  mod <- Phi <- as.list(1:d)
  eps <- matrix(0, nrow = n, ncol = d)
  for(i in 1:d){ ##fit univariate score models
    mod[[i]] <- ar(x[,i], order.max = p, demean = T, method = "ols", aic = F)
    Phi[[i]] <- mod[[i]]$x.intercept
    if(global_resids){
      eps[,i] <- eps_global[,i]
    } else {
      eps[,i] <- as.matrix(mod[[i]]$resid);
    }
    eps[,i][1:p] <- 1e-4 ##solve NA
    if(p==1) Phi[[i]] <- matrix(c(Phi[[i]],mod[[i]]$ar[,,1]),1,2 ) #cbind(Phi,  matrix( mod$ar, nrow=d, ncol=d))
    if(p>1){ 
      for (jj in 1:p){ #collect parameters into mat
        Phi[[i]] <- cbind(Phi[[i]],  mod[[i]]$ar[jj,,])
      } 
    } 
  }

  if(do_bootstrap == "regression"){ ## RUN BOOTSTRAP
    #bsmod <- bsPhi <- as.list(1:d)
    #N <- 250
    bsResid <- matrix(0, nrow = 2*G, ncol = d)
    for(i in 1:d){ ##obtain residuals from start and end G-segments
      #bsResid[,i] <- rbind( as.matrix(ar(x[1:G,i], order.max = p, demean = T, method = "ols", aic = F)$resid),
      #                    as.matrix(ar(x[(n-G+1):n,i], order.max = p, demean = T, method = "ols", aic = F)$resid)  ) ##model on start/end
      temp_resids <- as.matrix(ar(x[,i], order.max = p, demean = T, method = "ols", aic = F)$resid)
      bsResid[,i] <- temp_resids[c(1:G, (n-G+1):n)]
    }
    bsResid <- na.omit(bsResid)
    n_tilde <- nrow(bsResid)
    max_m <- bootstrap(x,p,G, Phi,bsResid , n_tilde, M, estim, varEstim)
    D_n <- quantile(max_m, 1-alpha)
  } else {
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (d*p+1)/2 * log(log(n/G)) - log(2/3 * gamma( (d*p+1)/2)) ##CORRECTED
  D_n <- (b+c_alpha)/a #threshold
  D_n <- max(D_n, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  }
  Reject <- FALSE
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector

  if(method == "Wald"){
    for (i in 1:d) {
      stat <- stat +  get_W_RCPP( as.matrix(x[,i]),p,G,estim)#W[,i] <- 
    }
  }
  if(method == "Score"){
  #if(global_resids) eps_global <-  (ar(x, order.max = p, demean = T, method = "ols", aic = F)$resid)
    stat <- get_T_RCPP(as.matrix(x), p, G, Phi = matrix(0), eps, PhiList = Phi, var_estim = varEstim, univariate = T )
  }
  
  cps <- c() #empty changepoint vector
  if(max(stat) > D_n){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(stat,D_n,G, nu=nu, criterion)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  } 
  
  ##Multiplier Bootstrap--------------------------
  if (do_bootstrap == "multiplier"){
    if( is.null(cps)) cps <- c(0)
    mbs <- multiplier_bootstrap(x, p, G, Phi, eps, cps, L = floor(n/4), M, estim, varEstim)
    D_n <- quantile(mbs, 1-alpha) ##overwrite threshold with bootstrap quantile 
    if(max(stat) > D_n){ #compare test stat with new threshold
      Reject <- TRUE
      cps <- get_cps(stat,D_n,G, nu=nu, criterion)
      if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
    } 
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

# mosum_univ(dp2_change, p=2, G=200)
# mosum_univ(dp2_change, p=2, G=200, method = "Score")
# mosum_univ(dp2_change, p=2, G=200, method = "Score", rm_cross_terms = T, do_bootstrap = "multiplier")
# 
# mosum_univ(dp2_change, p=2, G=200, method = "Score", rm_cross_terms = F, do_bootstrap = "multiplier")


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
# remove_cross_terms(dp2_change, 4, 3)

# 
# mpbootstrap <- function(x, p, cps, L, M){
#   n <- nrow(x)
#   K <- floor(n/L)
#   starts <- floor(L*cps/n); ends <- ceiling(L*cps/n) #isolate blocks containing cps
#   ## pre-process x
#   xbs <- x
#   for (ii in 1:length(cps)) {
#     indexset <- (K*starts[ii]):(K*ends[ii])
#     xbs[indexset,] <- colMeans(xbs[indexset,])
#   }
#   
#   Tm <- rep(0, M)
#   for (m in 1:M) {
#     perturb <- rep_len(rnorm(L), each = K, length.out = n)
#     xm <- perturb * xbs
#     Tm[m] <- get_T_RCPP()
#   }
#   return(xbs)
# }
# mpbootstrap(dp2_change, p=2, cps = c(1000,2000), L = 100, M =1)


# 
# bootstrap <- function(x_i, p, N) {
#   n <- length(x_i)
#   mod <- ar.ols(x_i[1:N,], aic = F, order.max = p)
#   e <- na.omit(mod$resid)
#   interval <- 1:nrow(e)
#   s <- sample(interval, size = n, replace = T)
#   bs <- as.matrix(e[s])
#   return(bs)
# }

# mbps_list <- list( matrix(c(0,0.5), 1,2 ),matrix(c(0,0.5), 1,2 ), matrix(c(0,0.5), 1,2 ))
##multiplier_bootstrap(dp2_change, 1, 200, mbps_list, dp2_eps, c(1), 100, 100, "DiagC", "Local")

sim_bootstrap <- function(x, p, G, method = "Score", estim = "DiagC", varEstim = "Local",  alpha = 0.05, criterion="eps", nu=.25,
                          rm_cross_terms =F, do_bootstrap = "regression", M = 1000, global_resids = F){
  model <- (ar.ols(x, F, p))
  model_list <- list()
  for (i in 1:p) {
    model_list[[i]] <-  (model$ar[i,,])
  }
  n <- nrow(x)
  d <- ncol(x)
  
  
  TM <- rep(0, M)
  r0 <- rep(0,d)
  Sigma <- model$var.pred #diag(1, d,d)
  for (m in 1:M) {
    xm <- VAR_sim(n, r0, Sigma,model_list, "normal", as.matrix(1), as.matrix(1)) #+ model$x.intercept
    TM[m] <- max(mosum_univ(xm, p, G, "Score")$stat)
  }
  return(TM)
}
sim_bootstrap(dp2_change, 2, 100, M = 20)
