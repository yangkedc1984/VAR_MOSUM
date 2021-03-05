library(devtools)
library(Matrix)
library(Rcpp)
library(mosumvar)

simnGdp <- function(iterations = 1000, n, G, p, test = "Score", estim = "DiagC", A1, criterion = "eta", nu=0.5){ 
  out <- vector(length = iterations)
  d <- nrow(A1)
  for (ii in 1:iterations) {
    
    #e1 <- matrix(rnorm(d* n, 0, .5),ncol=d) 
    x1 <-VAR_sim(n = n, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = list(A1), error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    
    ##null
    #sim_n41 <- as.matrix(rSim(A1,e1))#rbind(rSim_p(coeff_list1,e1),rSim_p(coeff_list2,e2),rSim_p(coeff_list3,e3)) #H0 or H1
    if(test =="Score")t_n41 <- test_Score_new(x=x1, p=p, G=G, Phi = NULL, eps = NULL, alpha = 0.05, estim=estim, var_estim = "Local", criterion, nu) 
    if(test =="Wald") t_n41 <- test_Wald_new(x=x1, p=p, G, alpha = 0.05, estim, criterion, nu)
    
    gc()
    out[ii] <- t_n41$Reject
  }
  return(out)
}

getMS <- function(rho, d){
  MS <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      MS[i,j] <- rho^(2 + 2*abs(i-j) )
    }
  }
  return(MS)
}


n <- 2000; d <- 3; p <- 1; Delta <- d*(d*p+1)
AnGdp <- getMS(0.5, d)
# 
# nGdpScore_3_d3 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps")
# nGdpScore_5_d3 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Score", A1 = AnGdp)
# nGdpScore_7_d3 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Score", A1 = AnGdp)
# 
# nGdpWald_3_d3 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Wald", A1 = AnGdp)
# nGdpWald_5_d3 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Wald", A1 = AnGdp)
# nGdpWald_7_d3 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Wald", A1 = AnGdp)
# save.image("nGdp.Rdata")
# 
# mean(nGdpScore_3_d3)
# mean(nGdpScore_5_d3)
# mean(nGdpScore_7_d3)
# mean(nGdpWald_3_d3)
# mean(nGdpWald_5_d3)
# mean(nGdpWald_7_d3)
# 
# 
# d <- 5; Delta <- d*(d*p+1); AnGdp <- getMS(0.5, d)
# nGdpScore_3_d5 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Score", A1 = AnGdp,  criterion = "eps", nu=0.25)
# nGdpScore_5_d5 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)
# nGdpScore_7_d5 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)
# 
# nGdpWald_3_d5 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Wald", A1 = AnGdp)
# nGdpWald_5_d5 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Wald", A1 = AnGdp)
# nGdpWald_7_d5 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Wald", A1 = AnGdp)
# save.image("nGdp.Rdata")
# 
# nGdpScore_9_d5 <- simnGdp(1000, n, G = 9 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)
# nGdpScore_11_d5 <- simnGdp(1000, n, G = 11 * Delta, p, test = "Score", A1 = AnGdp)
# nGdpScore_13_d5 <- simnGdp(1000, n, G = 13 * Delta, p, test = "Score", A1 = AnGdp)
# 
# 
# mean(nGdpScore_3_d5)
# mean(nGdpScore_5_d5)
# mean(nGdpScore_7_d5)
# mean(nGdpWald_3_d5)
# mean(nGdpWald_5_d5)
# mean(nGdpWald_7_d5)



d <- 7; Delta <- d*(d*p+1); AnGdp <- getMS(0.5, d)
nGdpScore_3_d7 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)
nGdpScore_5_d7 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)
nGdpScore_7_d7 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)
nGdpScore_9_d7 <- simnGdp(1000, n, G = 9 * Delta, p, test = "Score", A1 = AnGdp, criterion = "eps", nu=0.25)


nGdpWald_3_d7 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Wald", A1 = AnGdp, criterion = "eps", nu=0.25)
nGdpWald_5_d7 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Wald", A1 = AnGdp, criterion = "eps", nu=0.25)
nGdpWald_7_d7 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Wald", A1 = AnGdp, criterion = "eps", nu=0.25)
nGdpWald_9_d7 <- simnGdp(1000, n, G = 9 * Delta, p, test = "Wald", A1 = AnGdp, criterion = "eps", nu=0.25)

mean(nGdpScore_3_d7)
mean(nGdpScore_5_d7)
mean(nGdpScore_7_d7)
mean(nGdpScore_9_d7)

mean(nGdpWald_3_d7)
mean(nGdpWald_5_d7)
mean(nGdpWald_7_d7)
mean(nGdpWald_9_d7)
save.image("nGdp7.Rdata")



# 
# 
# 
# d <- 10; Delta <- d*(d*p+1); AnGdp <- getMS(0.5, d)
# nGdpScore_3_d10 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Score", A1 = AnGdp)
# nGdpScore_5_d10 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Score", A1 = AnGdp)
# nGdpScore_7_d10 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Score", A1 = AnGdp)
# 
# nGdpWald_3_d10 <- simnGdp(1000, n, G = 3 * Delta, p, test = "Wald", A1 = AnGdp)
# nGdpWald_5_d10 <- simnGdp(1000, n, G = 5 * Delta, p, test = "Wald", A1 = AnGdp)
# nGdpWald_7_d10 <- simnGdp(1000, n, G = 7 * Delta, p, test = "Wald", A1 = AnGdp) #dlog(d)?
# save.image("nGdp.Rdata")
# 
# 
# mean(nGdpScore_3_d10)
# mean(nGdpScore_5_d10)
# mean(nGdpScore_7_d10)
# mean(nGdpWald_3_d10)
# mean(nGdpWald_5_d10)
# mean(nGdpWald_7_d10)