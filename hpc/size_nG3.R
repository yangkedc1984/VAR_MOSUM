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

MS1d3 <- getMS(-0.6, 3)
# MS2d3 <- getMS(0.5, 3)
# MS3d3 <- getMS(-0.4, 3)
# 
MS1d5 <- getMS(-0.6, 5)
# MS2d5 <- getMS(0.5, 5)
# MS3d5 <- getMS(-0.4, 5)

MS1d7 <- getMS(-0.6, 7)
# MS2d7 <- getMS(0.5, 3)
# MS3d7 <- getMS(-0.4, 3)

MS1d10 <- getMS(-0.6, 10)
# MS2d10 <- getMS(0.5, 10)
# MS3d10 <- getMS(-0.4, 10)

MS1d20 <- getMS(-0.6, 20)


simSize <- function(iterations = 1000, n = 2000, G, p, test = "Score", estim = "DiagC", A1,A2=A1,A3=A1, truelag = 1, criterion = "eta", nu=0.5, option = F,
                    rm_cross_terms =T, global_resids = T){ #conservative nu
  coeff_list1 <-coeff_list2 <-coeff_list3 <- as.list(1:truelag) #empty AR param lists
  d <- nrow(A1)
  for (i in 1:truelag) {
    coeff_list1[[i]] <- A1 ^i#%^% i ##element-wise powers? division?
    #   coeff_list2[[i]] <- A2 ^i#%^% i
    #   coeff_list3[[i]] <- A3 ^i#%^% i
  }
  out <- vector("logical",iterations)
  d <- nrow(A1)
  for (ii in 1:iterations) {
    
    
    x1 <- mosumvar::VAR_sim(n = n, mu = rep(0,d), Sigma = diag(0.5,d,d), coeffs = coeff_list1, error_dist = "normal", P1 = matrix(1), Q1 =  matrix(1) )
    sim_n41 <-as.matrix(x1)
    
    
    if(option == "ar1"){
      t_n41 <- mosumvar::mosum_univ(sim_n41,p,G,method = test,  criterion = criterion, nu = nu,
                                    rm_cross_terms =rm_cross_terms, global_resids = global_resids ,do_bootstrap =F, M = 0)
      
    } else {
      t_n41 <-mosumvar::mosumvar(sim_n41, p, G, method = test)
    }
    #bound <- n/50
    #int1000 <- t_n41$cps[t_n41$cps <= 1060 & t_n41$cps >= 940]
    #int2000 <- t_n41$cps[t_n41$cps <= 2060 & t_n41$cps >= 1940]
    gc()
    out[ii] <- max(t_n41$mosum)
  }
  return(out)
}

# 
# #n = 1000
# n <- 1000
# 
#   G<- n/8
#   n1000d3G8 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n1000d5G8 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n1000d7G8 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n1000d10G8 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n1000d20G8 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")
#   
#   G<- n/4
#   n1000d3G4 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n1000d5G4 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n1000d7G4 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n1000d10G4 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n1000d20G4 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")

#n = 2000
# n <- 2000
# 
#   G<- n/16
#   n2000d3G16 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n2000d5G16 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n2000d7G16 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n2000d10G16 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n2000d20G16 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")
#   
#   G<- n/8
#   n2000d3G8 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n2000d5G8 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n2000d7G8 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n2000d10G8 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n2000d20G8 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")
#   
#   G<- n/4
#   n2000d3G4 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n2000d5G4 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n2000d7G4 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n2000d10G4 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n2000d20G4 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")

#n = 4000
# n <- 4000
# 
#   G<- n/16
#   n4000d3G16 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n4000d5G16 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n4000d7G16 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n4000d10G16 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n4000d20G16 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")
#   
#   G<- n/8
#   n4000d3G8 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n4000d5G8 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n4000d7G8 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n4000d10G8 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n4000d20G8 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")
#   
#   G<- n/4
#   n4000d3G4 <- simSize(1000, n, G, p=1, A1=MS1d3)
#   n4000d5G4 <- simSize(1000, n, G, p=1, A1=MS1d5)
#   n4000d7G4 <- simSize(1000, n, G, p=1, A1=MS1d7)
#   n4000d10G4 <- simSize(1000, n, G, p=1, A1=MS1d10, option = "ar1")
#   n4000d20G4 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")

# #n = 8000
n <- 8000

  G<- n/16
  # n8000d3G16 <- simSize(1000, n, G, p=1, A1=MS1d3)
  # n8000d5G16 <- simSize(1000, n, G, p=1, A1=MS1d5)
  n8000d7G16 <- simSize(1000, n, G, p=1, A1=MS1d7)
  #n8000d10G16 <- simSize(1000, n, G, p=1, A1=MS1d10)#, option = "ar1")
 # n8000d20G16 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")

  G<- n/8
  #n8000d3G8 <- simSize(1000, n, G, p=1, A1=MS1d3)
  #n8000d5G8 <- simSize(1000, n, G, p=1, A1=MS1d5)
  n8000d7G8 <- simSize(1000, n, G, p=1, A1=MS1d7)
  #n8000d10G8 <- simSize(1000, n, G, p=1, A1=MS1d10)#, option = "ar1")
  #n8000d20G8 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")

  G<- n/4
  #n8000d3G4 <- simSize(1000, n, G, p=1, A1=MS1d3)
#n8000d5G4 <- simSize(1000, n, G, p=1, A1=MS1d5)
  n8000d7G4 <- simSize(1000, n, G, p=1, A1=MS1d7)
 # n8000d10G4 <- simSize(1000, n, G, p=1, A1=MS1d10)#, option = "ar1")
 # n8000d20G4 <- simSize(1000, n, G, p=1, A1=MS1d20, option = "ar1")
  
save.image("simSize_nG3_7.Rdata")




