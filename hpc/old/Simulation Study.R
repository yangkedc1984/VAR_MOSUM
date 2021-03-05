A1 <- matrix( c(.7, -.1, -.1, -.1,
                -.1,.7,  -.1, -.1,
                -.1, -.1, .7, -.1,
                -.1, -.1, -.1,.7), nrow = 4, ncol = 4 ) 
A2 <- matrix( c(.6, .1, .1, .1,
                .1,.6,  .1, .1,
                .1, .1, .6, .1,
                .1, .1, .1,.6), nrow = 4, ncol = 4 ) 
A3 <- matrix( c(.5, .1, .1, -.1,
                .1,.5,  -.1, .1,
                .1, -.1, .5, .1,
                -.1, .1, .1,.5), nrow = 4, ncol = 4 ) 
Reject_count <- 0
Rej_vec <- rep(FALSE, 100)
#for (replicate in 1:100) {
testloop <- function(d, G, test = "Score", estim = "DiagC"){
  e1 <- matrix(rnorm(4 * 500, 0, .5),ncol=4)
  e2 <- matrix(rnorm(4 * 500, 0, .5),ncol=4)
  e3 <- matrix(rnorm(4 * 500, 0, .5),ncol=4)
  e4 <- matrix(rnorm(4 * 500, 0, .5),ncol=4)
  ##null
  sim_n41 <- rbind(rSim(A1,e1),rSim(A1,e2),rSim(A1,e3),rSim(A1,e4)) #H0 or H1
  plot.ts(sim_n41)
  m_n41 <- ar(sim_n41, order.max = 1, demean = T, method = "ols") #remove mean for E[Y]=0
  #var_change_0 <- as.matrix(t( t(var_change) - (change_model$x.mean))) #centre
  m_n41_a <- cbind(m_n41$x.intercept, matrix(m_n41$ar, nrow=4, ncol=4))
  m_n41_res <- m_n41$resid
  if(test =="Score") t_n41 <- test_Score(x=sim_n41, p=1, G, Phi = m_n41_a, eps = m_n41_res, alpha = 0.05, estim)
  if(test =="Wald") t_n41 <- test_Wald(x=sim_n41, p=1, G, alpha = 0.05, estim)
  int500 <- t_n41$cps[t_n41$cps <= 540 & t_n41$cps >= 460]
  int1000 <- t_n41$cps[t_n41$cps <= 1040 & t_n41$cps >= 960]
  int1500 <- t_n41$cps[t_n41$cps <= 1540 & t_n41$cps >= 1460]
  gc()
  return(c(t_n41$Reject, length(t_n41$cps), length(int500), length(int1000), length(int1500) ))
}
library(parallel)
Rej_list <- mcmapply(1:50, FUN=testloop, MoreArgs = list(G=100 , test="Wald", estim="DiagH"),  mc.cores = getOption("mc.cores", 4L))#mc.cores = getOption("mc.cores", 4L))
rowMeans(Rej_list); sd(Rej_list[2,])
gc() # garb 
