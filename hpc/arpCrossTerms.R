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


simMisspec <- function(iterations = 1000, G, p, test = "Score", estim = "DiagC", A1,A2,A3, truelag = 1, criterion = "eta", nu=0.5, option = "truncate", rm_cross_terms =F, global_resids = F){ #conservative nu
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
    sim_n41 <-as.matrix(rbind(x1,x2,x3))


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
    if(option == "ar1")   t_n41 <- mosum_univ(sim_n41,p,G,method = test,  criterion = criterion, nu = nu,rm_cross_terms =rm_cross_terms, global_resids = global_resids)

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

#MSWaldL1G200 <-  simMisspec(3, G=200,p=1, test = "Score", estim = "DiagC",MS1,MS2,MS3, truelag = 1, option = "ar1")

# ### no treatment
#
#
ar1ScoreL1G200CrossTerms <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL1G400CrossTerms <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL1G600CrossTerms <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

ar1ScoreL2G200CrossTerms <-  simMisspec(1000, G=200,p=2, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL2G400CrossTerms <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL2G600CrossTerms <-  simMisspec(1000, G=600,p=2, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

ar1ScoreL3G200CrossTerms <-  simMisspec(1000, G=200,p=3, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL3G400CrossTerms <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL3G600CrossTerms <-  simMisspec(1000, G=600,p=3, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

ar1ScoreL4G200CrossTerms <-  simMisspec(1000, G=200,p=4, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL4G400CrossTerms <-  simMisspec(1000, G=400,p=4, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL4G600CrossTerms <-  simMisspec(1000, G=600,p=4, test = "Score", estim = "DiagC", MS1,MS2,MS3, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

# ###
ar1ScoreL1G200NullCrossTerms <-  simMisspec(1000, G=200,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL1G400NullCrossTerms <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL1G600NullCrossTerms <-  simMisspec(1000, G=600,p=1, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

ar1ScoreL2G200NullCrossTerms <-  simMisspec(1000, G=200,p=2, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL2G400NullCrossTerms <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL2G600NullCrossTerms <-  simMisspec(1000, G=600,p=2, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

ar1ScoreL3G200NullCrossTerms <-  simMisspec(1000, G=200,p=3, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL3G400NullCrossTerms <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL3G600NullCrossTerms <-  simMisspec(1000, G=600,p=3, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

ar1ScoreL4G200NullCrossTerms <-  simMisspec(1000, G=200,p=4, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL4G400NullCrossTerms <-  simMisspec(1000, G=400,p=4, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)
ar1ScoreL4G600NullCrossTerms <-  simMisspec(1000, G=600,p=4, test = "Score", estim = "DiagC", MS1,MS1,MS1, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, rm_cross_terms = T)

save.image("arpCrossTerms.Rdata")


##########
report(ar1ScoreL1G200CrossTerms)
report(ar1ScoreL1G400CrossTerms)
report(ar1ScoreL1G600CrossTerms)

report(ar1ScoreL2G200CrossTerms)
report(ar1ScoreL2G400CrossTerms)
report(ar1ScoreL2G600CrossTerms)

report(ar1ScoreL3G200CrossTerms)
report(ar1ScoreL3G400CrossTerms)
report(ar1ScoreL3G600CrossTerms)

report(ar1ScoreL4G200CrossTerms)
report(ar1ScoreL4G400CrossTerms)
report(ar1ScoreL4G600CrossTerms)



mean(ar1ScoreL1G200NullCrossTerms[,1])
mean(ar1ScoreL1G400NullCrossTerms[,1])
mean(ar1ScoreL1G600NullCrossTerms[,1])

mean(ar1ScoreL2G200NullCrossTerms[,1])
mean(ar1ScoreL2G400NullCrossTerms[,1])
mean(ar1ScoreL2G600NullCrossTerms[,1])

mean(ar1ScoreL3G200NullCrossTerms[,1])
mean(ar1ScoreL3G400NullCrossTerms[,1])
mean(ar1ScoreL3G600NullCrossTerms[,1])

mean(ar1ScoreL4G200NullCrossTerms[,1])
mean(ar1ScoreL4G400NullCrossTerms[,1])
mean(ar1ScoreL4G600NullCrossTerms[,1])
