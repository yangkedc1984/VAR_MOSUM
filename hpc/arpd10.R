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

 MS1d3 <- diag(0.0, 3)# getMS(-0.6, 3)
 MS2d3 <- diag(0.0, 3)#getMS(0.5, 3)
 MS3d3 <- diag(0.0, 3)#getMS(-0.4, 3)
 
 MS1d3[row(MS1d3) == (col(MS1d3) - 1)] <- 0.6
 MS2d3 <- -MS1d3
 MS3d3 <- MS1d3
  # 
 
 
 MS1d5 <- getMS(-0.6, 5)
 MS2d5 <- getMS(0.5, 5)
 MS3d5 <- getMS(-0.4, 5)
 # MS1d5 <- diag(0.1, 5)# getMS(-0.6, 3)
 # MS2d5 <- diag(0.1, 5)#getMS(0.5, 3)
 # MS3d5 <- diag(0.1, 5)#getMS(-0.4, 3)
 
 MS1d5[row(MS1d5) == (col(MS1d5) - 1)] <- 0.9
 MS2d5 <- -MS1d5
 #MS2d5[row(MS2d5) == (col(MS2d5) - 1)] <- -0.9
 MS3d5 <- MS1d5
 
MS1d7 <- matrix(rnorm(49,0,0.1),7,7)# getMS(-0.6, 3)
MS2d7 <- matrix(rnorm(49,0,0.1),7,7)#getMS(0.5, 3)
MS3d7 <- matrix(rnorm(49,0,0.1),7,7)#getMS(-0.4, 3)

diag(MS1d7) <- 0.6
diag(MS2d7) <- 0.4
diag(MS3d7) <- 0.6
 
MS1d10 <- getMS(-0.6, 10)
MS2d10 <- getMS(0.5, 10)
MS3d10 <- getMS(-0.4, 10)


MS1d15 <- diag(0.1, 15)
MS2d15 <- diag(0.75, 15)
MS3d15 <- diag(-0.8, 15)


MS1d20 <- diag(0, 20)
MS1d20[row(MS1d20) == (col(MS1d20) - 1)] <- 0.6
MS2d20 <- - MS1d20
MS3d20 <- MS1d20

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
    x1 <- mosumvar::VAR_sim(n = 300, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list1, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    x2 <- mosumvar::VAR_sim(n = 300, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list2, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    x3 <- mosumvar::VAR_sim(n = 300, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list3, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
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
    if(option == "ar1")   t_n41 <- mosum_univ(sim_n41,p,G,method = test,  criterion = criterion, nu = nu,
                                              rm_cross_terms =rm_cross_terms, global_resids = global_resids ,do_bootstrap =do_bootstrap, M = M)

    # if(test=="BS") {
    #   tt <- MOSUMBS(x=sim_n41, p=lag, G)
    #   cps <- tt$cps
    #   if(is.null(cps)) cps <- c()
    #   t_n41 <- list(cps = cps, Reject = length(cps)>0)
    # }
    int1000 <- t_n41$cps[t_n41$cps <= 318 & t_n41$cps >= 282]
    int2000 <- t_n41$cps[t_n41$cps <= 618 & t_n41$cps >= 582]
    gc()
    out[ii,] <- c(t_n41$Reject, length(t_n41$cps), length(int1000), length(int2000) )
    #save(out, file ="arpd10.Rdata")
  }
  return(out)
}
##vary rho
diag15 <- diag(1,15)

null.05 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.05*diag15,0.05*diag15,0.05*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.10 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.1*diag15,0.1*diag15,0.1*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.15 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.15*diag15,0.15*diag15,0.15*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.20 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.20*diag15, 0.20*diag15,0.20*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.25 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.25*diag15,0.25*diag15,0.25*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 

alt.05 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.05*diag15,-0.05*diag15,0.05*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.10 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.1*diag15,-0.1*diag15,0.1*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.15 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.15*diag15,-0.15*diag15,0.15*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.20 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.20*diag15,-0.20*diag15,0.20*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.25 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", 0.25*diag15,-0.25*diag15,0.25*diag15, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 

mean(null.05[,1])
mean(null.10[,1])
mean(null.15[,1])
mean(null.20[,1])
mean(null.25[,1])

colMeans(alt.05)
colMeans(alt.10)
colMeans(alt.15)
colMeans(alt.20)
colMeans(alt.25)

sd(alt.05[,2])
sd(alt.10[,2])
sd(alt.15[,2])
sd(alt.20[,2])
sd(alt.25[,2])

## vary d
diag3 <- diag(.15,3)
diag5 <- diag(.15,5)
diag7 <- diag(.15,7)
diag10 <- diag(.15,10)

null.d3 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag3,diag3,diag3, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.d5 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag5,diag5,diag5, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.d7 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag7,diag7,diag7, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
null.d10 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag10,diag10,diag10, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 

alt.d3 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag3,-diag3,diag3, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.d5 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag5,-diag5,diag5, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.d7 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag7,-diag7,diag7, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 
alt.d10 <- simMisspec(1000, G=150,p=1, test = "Score", estim = "DiagC", diag10,-diag10,diag10, truelag = 1, option = "ar1",criterion = "eps", nu=0.2, global_resids = T, rm_cross_terms = T, do_bootstrap =F, M=300) 


