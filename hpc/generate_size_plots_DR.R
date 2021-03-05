# #load("~/Documents/Changepoints/hpc/simSize_nG.Rdata")
# c_alpha <- -log(log( (1-.05)^(-1/2))) 
# 
qgumDR <- function(x, n, G, d){
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (d+1)/2 * log(log(n/G)) - log(2/3 * gamma((d+1)/2))
  D_n <- (b+c_alpha)/a #threshold
  enforce <- sqrt(2*log(n)) + c_alpha/sqrt(2*log(n))
  D_n <- max(D_n, enforce)
  # if(D_n > enforce){
  #   p <- (x+b)/a#a*x-b
  #   return(1- exp(-2*exp(-p)) )
  # } else {
  #   return(1- exp(-2*exp(-x)))
  # }
  return(x/D_n)
}

library(ggplot2)


q1000DR <- as.data.frame( rbind( cbind(Gn = 1/8, d = 10, q= qgum(n1000d10G8, 1000, 1000/8, 10 )),
                cbind(Gn = 1/8, d = 15, q= qgum(n1000d15G8, 1000, 1000/8, 15 )),
                cbind(Gn = 1/8, d = 20, q= qgum(n1000d20G8, 1000, 1000/8, 20)),
                cbind(Gn = 1/4, d = 10, q= qgum(n1000d10G4, 1000, 1000/4, 10 )),
                cbind(Gn = 1/4, d = 15, q= qgum(n1000d15G4, 1000, 1000/4, 15 )),
                cbind(Gn = 1/4, d = 20, q= qgum(n1000d20G4, 1000, 1000/4, 20)) )
                )
q1000DR$Gn <- factor(q1000DR$Gn)
q1000DR$d <- factor(q1000DR$d)

par(mfrow = c(1,2))
p1000DR <- ggplot(data = q1000DR, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p1000DRplus <- p1000DR + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(title="n=1000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())


# df1000 <- data.frame(n = rep(1000, 8), G = c(rep(125,4),rep(250,4) ), size = n1000, d = rep(c(3,5,7,10),2) )
# df1000$Gn <- df1000$G/df1000$n 
# barplot( df1000$size, names.arg = df1000$d, xlab = "dimension", ylab = "size", main = "n=1000", col= 16*df1000$Gn )
# abline(h = 0.05, lty = 2, col = "red")


q2000DR <-as.data.frame( rbind( cbind(Gn = 1/16, d = 10, q= qgum(n2000d10G16, 2000, 2000/16, 10 )),
                                cbind(Gn = 1/16, d = 15, q= qgum(n2000d15G16, 2000, 2000/16, 15 )),
                                cbind(Gn = 1/16, d = 20, q= qgum(n2000d20G16, 2000, 2000/16, 20)),
                              cbind(Gn = 1/8, d = 10, q= qgum(n2000d10G8, 2000, 2000/8, 10 )),
                              cbind(Gn = 1/8, d = 15, q= qgum(n2000d15G8, 2000, 2000/8, 15 )),
                              cbind(Gn = 1/8, d = 20, q= qgum(n2000d20G8, 2000, 2000/8, 20)),
                              cbind(Gn = 1/4, d = 10, q= qgum(n2000d10G4, 2000, 2000/4, 10 )),
                              cbind(Gn = 1/4, d = 15, q= qgum(n2000d15G4, 2000, 2000/4, 15 )),
                              cbind(Gn = 1/4, d = 20, q= qgum(n2000d20G4, 2000, 2000/4, 20)) )
)
q2000DR$Gn <- factor(q2000DR$Gn)
q2000DR$d <- factor(q2000DR$d)

par(mfrow = c(1,2))
p2000DR <- ggplot(data = q2000DR, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p2000DRplus <- p2000DR + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(title="n=2000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())


##
#load("~/Documents/Changepoints/hpc/simSize_nG2.Rdata")
q4000DR <-as.data.frame( rbind( cbind(Gn = 1/16, d = 10, q= qgum(n4000d10G16, 4000, 4000/16, 10 )),
                                cbind(Gn = 1/16, d = 15, q= qgum(n4000d15G16, 4000, 4000/16, 15 )),
                                cbind(Gn = 1/16, d = 20, q= qgum(n4000d20G16, 4000, 4000/16, 20)),
                                cbind(Gn = 1/8, d = 10, q= qgum(n4000d10G8, 4000, 4000/8, 10 )),
                                cbind(Gn = 1/8, d = 15, q= qgum(n4000d15G8, 4000, 4000/8, 15 )),
                                cbind(Gn = 1/8, d = 20, q= qgum(n4000d20G8, 4000, 4000/8, 20)),
                                cbind(Gn = 1/4, d = 10, q= qgum(n4000d10G4, 4000, 4000/4, 10 )),
                                cbind(Gn = 1/4, d = 15, q= qgum(n4000d15G4, 4000, 4000/4, 15 )),
                                cbind(Gn = 1/4, d = 20, q= qgum(n4000d20G4, 4000, 4000/4, 20))
) )
q4000DR$Gn <- factor(q4000DR$Gn)
q4000DR$d <- factor(q4000DR$d)

par(mfrow = c(1,2))
p4000DR <- ggplot(data = q4000DR, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p4000DRplus <- p4000DR + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(title="n=4000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

##
q8000DR <-as.data.frame( rbind( cbind(Gn = 1/16, d = 10, q= qgum(n8000d10G16, 8000, 8000/16, 10 )),
                                cbind(Gn = 1/16, d = 15, q= qgum(n8000d15G16, 8000, 8000/16, 15 )),
                                cbind(Gn = 1/16, d = 20, q= qgum(n8000d20G16, 8000, 8000/16, 20)),
                                cbind(Gn = 1/8, d = 10, q= qgum(n8000d10G8, 8000, 8000/8, 10 )),
                                cbind(Gn = 1/8, d = 15, q= qgum(n8000d15G8, 8000, 8000/8, 15 )),
                                cbind(Gn = 1/8, d = 20, q= qgum(n8000d20G8, 8000, 8000/8, 20)),
                                cbind(Gn = 1/4, d = 10, q= qgum(n8000d10G4, 8000, 8000/4, 10 )),
                                cbind(Gn = 1/4, d = 15, q= qgum(n8000d15G4, 8000, 8000/4, 15 )),
                                cbind(Gn = 1/4, d = 20, q= qgum(n8000d20G4, 8000, 8000/4, 20)) )
)
q8000DR$Gn <- factor(q8000DR$Gn)
q8000DR$d <- factor(q8000DR$d)

par(mfrow = c(1,2))
p8000DR <- ggplot(data = q8000DR, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p8000DRplus <- p8000DR + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(title="n=8000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

library(gridExtra)
grid.arrange(p1000DRplus, p2000DRplus, p4000DRplus, p8000DRplus)
