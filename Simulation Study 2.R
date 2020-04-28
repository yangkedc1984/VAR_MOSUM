A1_r1 <- matrix( c(.5, .1, .1, 
                 .1, .5,  .1, 
                 .1, .1, .5), nrow = 3, ncol = 3 ) 
A2_r1 <- matrix( c(-.2, .2, .2, 
                   .2, -.2,  .2, 
                   .2, .2, -.2), nrow = 3, ncol = 3 ) 

A1_r2 <- A2_r1
A2_r2 <- A1_r1

A1_r3 <- matrix( c(.4, .4, .0, 
                   .4, .4,  .4, 
                   .0, .4, .4), nrow = 3, ncol = 3 ) 
A2_r3 <- matrix( c(-.2, .0, .0, 
                   .0, -.2,  .0, 
                   .0, .0, -.2), nrow = 3, ncol = 3 ) 

  
Reject_count <- 0
Rej_vec <- rep(FALSE, 100)
#for (replicate in 1:100) {
testloop <- function(d, G, test = "Score", estim = "DiagC"){
  e1 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  e2 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  e3 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  e4 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  ##null
  sim_n41 <- rbind(rSim_p2(A1_r1,A2_r1,e1),rSim_p2(A1_r1,A2_r1,e2),rSim_p2(A1_r1,A2_r1,e3),rSim_p2(A1_r1,A2_r1,e4)) #H0 or H1
  if(test =="Score"){ 

    #plot.ts(sim_n41)
    m_n41 <- ar(sim_n41, order.max = 2, aic = F, demean = T, method = "ols") #remove mean for E[Y]=0
    #var_change_0 <- as.matrix(t( t(var_change) - (change_model$x.mean))) #centre
    m_n41_a <- cbind(m_n41$x.intercept, matrix(m_n41$ar[1,,], nrow=3, ncol=3),matrix(m_n41$ar[2,,], nrow=3, ncol=3))
    m_n41_res <- m_n41$resid
    t_n41 <- test_Score(x=sim_n41, p=2, G, Phi = m_n41_a, eps = m_n41_res, alpha = 0.05, estim)
  }
  if(test =="Wald") t_n41 <- test_Wald(x=sim_n41, p=2, G, alpha = 0.05, estim)
  int500 <- t_n41$cps[t_n41$cps <= 540 & t_n41$cps >= 460]
  int1000 <- t_n41$cps[t_n41$cps <= 1040 & t_n41$cps >= 960]
  int1500 <- t_n41$cps[t_n41$cps <= 1540 & t_n41$cps >= 1460]
  gc()
  return(c(t_n41$Reject, length(t_n41$cps), length(int500), length(int1000), length(int1500) ))
}

testloop_alt <- function(d, G, test = "Score", estim = "DiagC"){
  e1 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  e2 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  e3 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  e4 <- matrix(rnorm(3 * 500, 0, .5),ncol=3)
  ##alt
  sim_n41 <- rbind(rSim_p2(A1_r1,A2_r1,e1),rSim_p2(A1_r2,A2_r2,e2),rSim_p2(A1_r3,A2_r3,e3),rSim_p2(A1_r1,A2_r1,e4)) #H0 or H1
  if(test =="Score"){ 
    #plot.ts(sim_n41)
    m_n41 <- ar(sim_n41, order.max = 2, aic = F, demean = T, method = "ols") #remove mean for E[Y]=0
    #var_change_0 <- as.matrix(t( t(var_change) - (change_model$x.mean))) #centre
    m_n41_a <- cbind(m_n41$x.intercept, matrix(m_n41$ar[1,,], nrow=3, ncol=3),matrix(m_n41$ar[2,,], nrow=3, ncol=3))
    m_n41_res <- m_n41$resid
    t_n41 <- test_Score(x=sim_n41, p=2, G, Phi = m_n41_a, eps = m_n41_res, alpha = 0.05, estim)
  }
  if(test =="Wald") t_n41 <- test_Wald(x=sim_n41, p=2, G, alpha = 0.05, estim)
  int500 <- t_n41$cps[t_n41$cps <= 540 & t_n41$cps >= 460]
  int1000 <- t_n41$cps[t_n41$cps <= 1040 & t_n41$cps >= 960]
  int1500 <- t_n41$cps[t_n41$cps <= 1540 & t_n41$cps >= 1460]
  gc()
  return(c(t_n41$Reject, length(t_n41$cps), length(int500), length(int1000), length(int1500) ))
}
library(parallel)
Rej_list <- mcmapply(1:100, FUN=testloop, MoreArgs = list(G=200 , test="Score", estim="DiagH"),  mc.cores = getOption("mc.cores", 4L))#mc.cores = getOption("mc.cores", 4L))
rowMeans(Rej_list); sd(Rej_list[2,])

DC100 <- mcmapply(1:100, FUN=testloop, MoreArgs = list(G=100 , test="Wald", estim="DiagC"),  mc.cores = getOption("mc.cores", 4L))
DC150 <- mcmapply(1:100, FUN=testloop, MoreArgs = list(G=150 , test="Wald", estim="DiagC"),  mc.cores = getOption("mc.cores", 4L))
DC200 <- mcmapply(1:100, FUN=testloop, MoreArgs = list(G=200 , test="Wald", estim="DiagC"),  mc.cores = getOption("mc.cores", 4L))
DC100A <- mcmapply(1:100, FUN=testloop_alt, MoreArgs = list(G=100 , test="Wald", estim="DiagC"),  mc.cores = getOption("mc.cores", 4L))
DC150A <- mcmapply(1:100, FUN=testloop_alt, MoreArgs = list(G=150 , test="Wald", estim="DiagC"),  mc.cores = getOption("mc.cores", 4L))
DC200A <- mcmapply(1:100, FUN=testloop_alt, MoreArgs = list(G=200 , test="Wald", estim="DiagC"),  mc.cores = getOption("mc.cores", 4L))
gc() # garb 
