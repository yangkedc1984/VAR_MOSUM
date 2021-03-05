library(devtools)
library(Matrix)
library(Rcpp)
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

d<- 5
MS1 <- getMS(-0.6, d)
MS2 <- getMS(0.5, d)
MS3 <- getMS(-0.4, d)
P1 <- getMS(-0.5, d)
Q1 <- getMS(-0.5, d)



simNonnormal <- function(iterations = 1000, G=400, p, test = "Score", estim = "DiagC", error_dist, A1,A2,A3, criterion = "eta", nu=0.5, df=4){ #conservative nu
  coeff_list1 <-coeff_list2 <-coeff_list3 <- as.list(1:p) #empty AR param lists
  d <- nrow(A1)
  for (i in 1:p) {
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
    x1 <- VAR_sim(n = 1000, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list1, error_dist = error_dist, P1 =P1, Q1 = Q1, df )
    x2 <- VAR_sim(n = 1000, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list2, error_dist = error_dist, P1 =P1, Q1 = Q1, df)
    x3 <- VAR_sim(n = 1000, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list3, error_dist = error_dist, P1 =P1, Q1 = Q1, df )
    sim_n41 <-rbind(x1,x2,x3)
    
    
    if(test =="Score"){
      m_n41 <- ar(sim_n41, order.max = p, aic = F, demean = T, method = "ols")
      m_n41_a <- m_n41$x.intercept#cbind(m_n41$x.intercept, as.matrix(m_n41$ar[1,,]))
      for (i in 1:p) {
        m_n41_a <- cbind(m_n41_a, as.matrix(m_n41$ar[i,,]))
      }
      m_n41_res <- m_n41$resid; m_n41_res[1:d,] <- 0.0001
      t_n41 <- test_Score_new(x=sim_n41, p=p, G, Phi = as.matrix(m_n41_a), eps = as.matrix(m_n41_res), alpha = 0.05, estim, var_estim = "Local", criterion, nu)
    }  
    if(test =="Wald") t_n41 <- test_Wald_new(x=sim_n41, p=p, G, alpha = 0.05, estim, criterion, nu)
    
    
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

# simNonnormal(1,400,1,"Score","DiagC","normal",MS1,MS2,MS3)
# simNonnormal(1,400,1,"Score","DiagC","t",MS1,MS2,MS3, df=4)
# simNonnormal(1,400,1,"Score","DiagC","garch",MS1,MS2,MS3)

## normal
# normalScore <- simNonnormal(1000,400,1,"Score","DiagC","normal",MS1,MS2,MS3)
# normalWald <- simNonnormal(1000,400,1,"Wald","DiagC","normal",MS1,MS2,MS3)
# 
 normalScoreNull <- simNonnormal(1000,400,1,"Score","DiagC","normal",MS1,MS1,MS1, criterion = "eps", nu=0.25)
# normalWaldNull <- simNonnormal(1000,400,1,"Wald","DiagC","normal",MS1,MS1,MS1)

## t 4
#t4Score <- simNonnormal(1000,400,1,"Score","DiagC","t",MS1,MS2,MS3, df=4) ##RERUN
# t4Wald <- simNonnormal(1000,400,1,"Wald","DiagC","t",MS1,MS2,MS3, df=4)
# 
 t4ScoreNull <- simNonnormal(1000,400,1,"Score","DiagC","t",MS1,MS1,MS1, df=4, criterion = "eps", nu=0.25)
# t4WaldNull <- simNonnormal(1000,400,1,"Wald","DiagC","t",MS1,MS1,MS1, df=4)

## t 8
#t8Score <- simNonnormal(1000,400,1,"Score","DiagC","t",MS1,MS2,MS3, df=8) ##assign
# t8Wald <- simNonnormal(1000,400,1,"Wald","DiagC","t",MS1,MS2,MS3, df=8)
# 
 t8ScoreNull <- simNonnormal(1000,400,1,"Score","DiagC","t",MS1,MS1,MS1, df=8, criterion = "eps", nu=0.25)
# t8WaldNull <- simNonnormal(1000,400,1,"Wald","DiagC","t",MS1,MS1,MS1, df=8)


## BEKK garch
#garchScore <- simNonnormal(1000,400,1,"Score","DiagC","garch",MS1,MS2,MS3) ##RERUN
#garchWald <- simNonnormal(1000,400,1,"Wald","DiagC","garch",MS1,MS2,MS3) ##RERUN

 garchScoreNull <- simNonnormal(1000,400,1,"Score","DiagC","garch",MS1,MS1,MS1, criterion = "eps", nu=0.25) ##assign
 #garchWaldNull <- simNonnormal(1000,400,1,"Wald","DiagC","garch",MS1,MS1,MS1) ##assign

save.image("nonnormalEps.Rdata")
# 
# report(normalScore)
# report(normalWald)
# mean(normalScoreNull[,1])
# mean(normalWaldNull[,1])
# 
# report(t4Score)
# report(t4Wald)
# mean(t4ScoreNull[,1])
# mean(t4WaldNull[,1])
# 
# 
# report(t8Score)
# report(t8Wald)
# mean(t8ScoreNull[,1])
# mean(t8WaldNull[,1])
# 
# report(garchScore)
# report(garchWald)
# mean(garchScoreNull[,1])
# mean(garchWaldNull[,1])

