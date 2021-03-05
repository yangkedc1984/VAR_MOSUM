 library(devtools)
 library(Matrix)
 library(Rcpp)
 devtools::install_github("Dom-Owens-UoB/VAR_MOSUM/mosumvar")

 library(mosumvar)

getMS <- function(rho, d){
  MS <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      MS[i,j] <- rho^(2 + 2*abs(i-j) )
    }
  }
  return(MS)
}

# MS1d3 <- getMS(-0.6, 3)
# MS2d3 <- getMS(0.5, 3)
# MS3d3 <- getMS(-0.4, 3)
# 
# MS1d5 <- getMS(-0.6, 5)
# MS2d5 <- getMS(0.5, 5)
# MS3d5 <- getMS(-0.4, 5)

# MS1d7 <- getMS(-0.6, 3)
# MS2d7 <- getMS(0.5, 3)
# MS3d7 <- getMS(-0.4, 3)

MS1d20 <- getMS(-0.6, 20)
MS2d20 <- getMS(0.5, 20)
MS3d20 <- getMS(-0.4, 20)



simMisspec <- function(iterations = 1000, G, p, test = "Score", estim = "DiagC", A1,A2,A3, truelag = 1, criterion = "eta", nu=0.5, option = "truncate",
                       rm_cross_terms =F, global_resids = F, do_bootstrap = F, M = 1000){ #conservative nu
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
    if(option == "ar1")   t_n41 <- mosumvar::mosum_univ(sim_n41,p,G,method = test,  criterion = criterion, nu = nu,
                                              rm_cross_terms =rm_cross_terms, global_resids = global_resids ,do_bootstrap =do_bootstrap, M = M)

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
    #save(out, file ="arpd10Null.Rdata")
  }
  return(out)
}

#MSWaldL1G200 <-  simMisspec(3, G=200,p=1, test = "Score", estim = "DiagC",MS1,MS2,MS3, truelag = 1, option = "ar1")

# ### no treatment
#
# #
# d3Boot <-  simMisspec(500, G=400,p=1, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3NullBoot <-  simMisspec(500, G=400,p=1, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3Bootp2 <-  simMisspec(500, G=400,p=2, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 2, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3NullBootp2 <-  simMisspec(500, G=400,p=2, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3Bootp3 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 3, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3NullBootp3 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3Bootp4 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 4, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d3NullBootp4 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# 
#
#d5Boot <-  simMisspec(200, G=400,p=1, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
#d5NullBoot <-  simMisspec(200, G=400,p=1, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
#d5Bootp2 <-  simMisspec(500, G=400,p=2, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 2, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
#d5NullBootp2 <-  simMisspec(500, G=400,p=2, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d5Bootp3 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 3, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d5NullBootp3 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d5Bootp4 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 4, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d5NullBootp4 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)

# 
# d7Boot <-  simMisspec(500, G=400,p=1, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7NullBoot <-  simMisspec(500, G=400,p=1, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7Bootp2 <-  simMisspec(500, G=400,p=2, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 2, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7NullBootp2 <-  simMisspec(500, G=400,p=2, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7Bootp3 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 3, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7NullBootp3 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7Bootp4 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 4, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)
# d7NullBootp4 <-  simMisspec(500, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=500)

## no boot
# 
# d3no <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3Nullno<-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3nop2 <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 2, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3Nullnop2 <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3nop3 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 3, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3Nullnop3 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3nop4 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS2d3,MS3d3, truelag = 4, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d3Nullnop4 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d3,MS1d3,MS1d3, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# 
# 
# d5no <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5Nullno <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5nop2 <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 2, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5Nullnop2 <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5nop3 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 3, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5Nullnop3 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5nop4 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS2d5,MS3d5, truelag = 4, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d5Nullnop4 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d5,MS1d5,MS1d5, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# 
# 
# d7no <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7Nullno <-  simMisspec(1000, G=400,p=1, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 1, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7nop2 <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 2, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7Nullnop2 <-  simMisspec(1000, G=400,p=2, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 2, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7nop3 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 3, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7Nullnop3 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 3, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7nop4 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS2d7,MS3d7, truelag = 4, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)
# d7Nullnop4 <-  simMisspec(1000, G=400,p=3, test = "Score", estim = "DiagC", MS1d7,MS1d7,MS1d7, truelag = 4, option = "ar1", criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = F, M=500)





# d10BS <- simMisspec(1000, G=500,p=1, test = "Score", estim = "DiagC", MS1d10,MS2d10,MS3d10, truelag =1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=300)
# d10BS_Null <- simMisspec(1000, G=500,p=1, test = "Score", estim = "DiagC", MS1d10,MS1d10,MS1d10, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=300)

# d15BS <- simMisspec(1000, G=500,p=1, test = "Score", estim = "DiagC", MS1d15,MS2d15,MS3d15, truelag =1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=300)
# d15BS_Null <- simMisspec(1000, G=500,p=1, test = "Score", estim = "DiagC", MS1d15,MS1d15,MS1d15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=300)

#d20BS <- simMisspec(1000, G=500,p=1, test = "Score", estim = "DiagC", MS1d20,MS2d20,MS3d20, truelag =1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=300)
d20BS_Null <- simMisspec(1000, G=500,p=1, test = "Score", estim = "DiagC", MS1d20,MS1d20,MS1d20, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap = "multiplier", M=300)

save.image("arpd20BSNull.Rdata")

# 
 # report(d5Boot)
 # 
 # mean(d5NullBoot[,1])
# mean(ar1ScoreL2G400NullBoot)
# mean(ar1ScoreL3G400NullBoot)
# mean(ar1ScoreL4G400NullBoot)




