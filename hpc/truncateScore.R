 library(devtools)
 library(Matrix)
 library(Rcpp)
 library(mosumvar)
 #devtools::install_github("Dom-Owens-UoB/VAR_MOSUM/mosumvar")


getMS <- function(rho, d){
  MS <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      MS[i,j] <- rho^(2 + 2*abs(i-j) )
    }
  }
  return(MS)
}

MS1 <- getMS(-0.6, 10)
MS2 <- getMS(0.5, 10)
MS3 <- getMS(-0.4, 10)


simMisspec <- function(iterations = 1000, G, p, test = "Score", estim = "DiagC", A1,A2,A3, truelag = 1, criterion = "eta", nu=0.5, option = "truncate", rm_cross_terms =F){ #conservative nu
  coeff_list1 <-coeff_list2 <-coeff_list3 <- as.list(1:truelag) #empty AR param lists
  d <- nrow(A1)
  for (i in 1:truelag) {
    coeff_list1[[i]] <- A1 ^i#%^% i ##element-wise powers? division?
    coeff_list2[[i]] <- A2 ^i#%^% i
    coeff_list3[[i]] <- A3 ^i#%^% i
  }
  out <- matrix(NA,iterations,4)
  d <- nrow(A1)
  for (ii in 1:iterations) {

    # e1 <- matrix(rnorm(d* 1000, 0, .5),ncol=d) #n = 3000 here
    # e2 <- matrix(rnorm(d * 1000, 0, .5),ncol=d)
    # e3 <- matrix(rnorm(d * 1000, 0, .5),ncol=d)
    # ##null
    # sim_n41 <- rbind(rSim_p(coeff_list1,e1),rSim_p(coeff_list2,e2),rSim_p(coeff_list3,e3)) #H0 or H1
    x1 <- mosumvar::VAR_sim(n = 1000, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list1, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    x2 <- mosumvar::VAR_sim(n = 1000, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list2, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    x3 <- mosumvar::VAR_sim(n = 1000, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list3, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    sim_n41 <-rbind(x1,x2,x3)


    if(option == "truncate"){
      if(test == "Score"){
        m_n41 <- ar(sim_n41, order.max = p, aic = F, demean = T, method = "ols")
        m_n41_a <- m_n41$x.intercept#cbind(m_n41$x.intercept, as.matrix(m_n41$ar[1,,]))
        for (i in 1:p) {
          m_n41_a <- cbind(m_n41_a, as.matrix(m_n41$ar[i,,]))
        }
        m_n41_res <- m_n41$resid; m_n41_res[1:d,] <- 0.0001
        t_n41 <- test_Score_new(x=sim_n41, p=p, G, Phi = as.matrix(m_n41_a), eps = as.matrix(m_n41_res), alpha = 0.05, estim, var_estim = "Local", criterion, nu)
      }
      if(test =="Wald") t_n41 <- test_Wald_new(x=sim_n41, p=p, G, alpha = 0.05, estim, criterion, nu)
    }
    if(option == "ar1")   t_n41 <- mosum_univ(sim_n41,p,G,method = test, criterion = criterion, nu = nu,rm_cross_terms =rm_cross_terms)

    # if(test=="BS") {
    #   tt <- MOSUMBS(x=sim_n41, p=lag, G)
    #   cps <- tt$cps
    #   if(is.null(cps)) cps <- c()
    #   t_n41 <- list(cps = cps, Reject = length(cps)>0)
    # }
    int1000 <- t_n41$cps[t_n41$cps <= 1060 & t_n41$cps >= 940]
    int2000 <- t_n41$cps[t_n41$cps <= 2060 & t_n41$cps >= 1940]
    gc()
    out[ii,] <- c(t_n41$Reject, length(t_n41$cps), length(int1000), length(int2000) )
  }
  return(out)
}

#MSWaldL1G200 <-  simMisspec(3, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 1)


# ### Wald
# 
# MSWaldL1G200 <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 1)
# MSWaldL1G400 <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 1)
# MSWaldL1G600 <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 1)
# 
# MSWaldL2G200 <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 2)
# MSWaldL2G400 <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 2)
# MSWaldL2G600 <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 2)
# 
# MSWaldL3G200 <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 3)
# MSWaldL3G400 <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 3)
# MSWaldL3G600 <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 3)
# 
# MSWaldL4G200 <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 4)
# MSWaldL4G400 <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 4)
# MSWaldL4G600 <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS2,MS3, truelag = 4)
# 
# save.image("truncateWald.Rdata")
# 
# 
# # ###
# MSWaldL1G200Null <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 1)
# MSWaldL1G400Null <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 1)
# MSWaldL1G600Null <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 1)
# 
# MSWaldL2G200Null <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 2)
# MSWaldL2G400Null <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 2)
# MSWaldL2G600Null <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 2)
# 
# MSWaldL3G200Null <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 3)
# MSWaldL3G400Null <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 3)
# MSWaldL3G600Null <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 3)
# 
# MSWaldL4G200Null <-  simMisspec(1000, G=200,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 4)
# MSWaldL4G400Null <-  simMisspec(1000, G=400,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 4)
# MSWaldL4G600Null <-  simMisspec(1000, G=600,p=1, test = "Wald", estim = "DiagC", MS1,MS1,MS1, truelag = 4)
# 
# save.image("truncateWaldNull.Rdata")
# 
# #
# # ### Score
#
#
MSScoreL1G200 <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 1)
MSScoreL1G400 <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 1)
MSScoreL1G600 <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 1)

MSScoreL2G200 <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 2)
MSScoreL2G400 <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 2)
MSScoreL2G600 <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 2)

MSScoreL3G200 <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 3)
MSScoreL3G400 <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 3)
MSScoreL3G600 <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 3)

MSScoreL4G200 <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 4)
MSScoreL4G400 <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 4)
MSScoreL4G600 <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 4)

save.image("truncateScore.Rdata")
# 
# # ###
# MSScoreL1G200Null <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 1)
# MSScoreL1G400Null <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 1)
# MSScoreL1G600Null <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 1)
# 
# MSScoreL2G200Null <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 2)
# MSScoreL2G400Null <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 2)
# MSScoreL2G600Null <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 2)
# 
# MSScoreL3G200Null <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 3)
# MSScoreL3G400Null <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 3)
# MSScoreL3G600Null <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 3)
# 
# MSScoreL4G200Null <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 4)
# MSScoreL4G400Null <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 4)
# MSScoreL4G600Null <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 4)
# 
# save.image("truncateScoreNull.Rdata")

#save.image("truncate.Rdata")
