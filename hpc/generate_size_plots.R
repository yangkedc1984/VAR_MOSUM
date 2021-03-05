#load("~/Documents/Changepoints/hpc/simSize_nG.Rdata")
c_alpha <- -log(log( (1-.05)^(-1/2))) 

qgum <- function(x, n, G, d){
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + d*(d+1)/2 * log(log(n/G)) - log(2/3 * gamma(d*(d+1)/2))
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


boxplot(qgum(n1000d3G8, 1000, 1000/8, 3) )
boxplot( qgum(n1000d5G4, 1000, 1000/4, 5) )

q1000 <- as.data.frame( rbind( cbind(Gn = 1/8, d = 3, q= qgum(n1000d3G8, 1000, 1000/8, 3 )),
                cbind(Gn = 1/8, d = 5, q= qgum(n1000d5G8, 1000, 1000/8, 5 )),
                cbind(Gn = 1/8, d = 7, q= qgum(n1000d7G8, 1000, 1000/8, 7)),
                cbind(Gn = 1/8, d = 10, q= qgum(n1000d10G8, 1000, 1000/8, 10)),
                cbind(Gn = 1/4, d = 3, q=qgum(n1000d3G4, 1000, 1000/4, 3)),
                cbind(Gn = 1/4, d = 5, q=qgum(n1000d5G4, 1000, 1000/4, 5)),
                cbind(Gn = 1/4, d = 7, q=qgum(n1000d7G4, 1000, 1000/4, 7)),
                cbind(Gn = 1/4, d = 10, q=qgum(n1000d10G4, 1000, 1000/4, 10)) ) 
                )
q1000$Gn <- factor(q1000$Gn)
q1000$d <- factor(q1000$d)

par(mfrow = c(1,2))
p1000 <- ggplot(data = q1000, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p1000plus <- p1000 + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(title="n=1000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
save.image("~hpc/q1000.pdf")

# df1000 <- data.frame(n = rep(1000, 8), G = c(rep(125,4),rep(250,4) ), size = n1000, d = rep(c(3,5,7,10),2) )
# df1000$Gn <- df1000$G/df1000$n 
# barplot( df1000$size, names.arg = df1000$d, xlab = "dimension", ylab = "size", main = "n=1000", col= 16*df1000$Gn )
# abline(h = 0.05, lty = 2, col = "red")


q2000 <- as.data.frame( rbind(cbind(Gn = 1/16, d = 3, q=qgum(n2000d3G16, 2000, 2000/16, 3)),
                              cbind(Gn = 1/16, d = 5, q=qgum(n2000d5G16, 2000, 2000/16, 5)),
                              cbind(Gn = 1/16, d = 7, q=qgum(n2000d7G16, 2000, 2000/16, 7)),
                              cbind(Gn = 1/16, d = 10, q=qgum(n2000d10G16, 2000, 2000/16, 10)), 
                                cbind(Gn = 1/8, d = 3, q= qgum(n2000d3G8, 2000, 2000/8, 3 )),
                               cbind(Gn = 1/8, d = 5, q= qgum(n2000d5G8, 2000, 2000/8, 5 )),
                               cbind(Gn = 1/8, d = 7, q= qgum(n2000d7G8, 2000, 2000/8, 7)),
                               cbind(Gn = 1/8, d = 10, q= qgum(n2000d10G8, 2000, 2000/8, 10)),
                               cbind(Gn = 1/4, d = 3, q=qgum(n2000d3G4, 2000, 2000/4, 3)),
                               cbind(Gn = 1/4, d = 5, q=qgum(n2000d5G4, 2000, 2000/4, 5)),
                               cbind(Gn = 1/4, d = 7, q=qgum(n2000d7G4, 2000, 2000/4, 7)),
                               cbind(Gn = 1/4, d = 10, q=qgum(n2000d10G4, 2000, 2000/4, 10)) ) 
                        )
q2000$Gn <- factor(q2000$Gn)
q2000$d <- factor(q2000$d)

par(mfrow = c(1,3))
p2000 <- ggplot(data = q2000, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p2000plus <- p2000 + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red")+
  labs(title="n=2000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())


##
#load("~/Documents/Changepoints/hpc/simSize_nG2.Rdata")

q4000 <- as.data.frame( rbind(cbind(Gn = 1/16, d = 3, q=qgum(n4000d3G16, 4000, 4000/16, 3)),
                              cbind(Gn = 1/16, d = 5, q=qgum(n4000d5G16, 4000, 4000/16, 5)),
                              cbind(Gn = 1/16, d = 7, q=qgum(n4000d7G16_correct, 4000, 4000/16, 7)),
                              cbind(Gn = 1/16, d = 10, q=qgum(n4000d10G16, 4000, 4000/16, 10)), 
                              cbind(Gn = 1/8, d = 3, q= qgum(n4000d3G8, 4000, 4000/8, 3 )),
                              cbind(Gn = 1/8, d = 5, q= qgum(n4000d5G8, 4000, 4000/8, 5 )),
                              cbind(Gn = 1/8, d = 7, q= qgum(n4000d7G8_correct, 4000, 4000/8, 7)),
                              cbind(Gn = 1/8, d = 10, q= qgum(n4000d10G8_correct, 4000, 4000/8, 10)), ##
                              cbind(Gn = 1/4, d = 3, q=qgum(n4000d3G4, 4000, 4000/4, 3)),
                              cbind(Gn = 1/4, d = 5, q=qgum(n4000d5G4, 4000, 4000/4, 5)),
                              cbind(Gn = 1/4, d = 7, q=qgum(n4000d7G4_correct, 4000, 4000/4, 7)),
                              cbind(Gn = 1/4, d = 10, q=qgum(n4000d10G4, 4000, 4000/4, 10)) ) 
)
q4000$Gn <- factor(q4000$Gn)
q4000$d <- factor(q4000$d)

par(mfrow = c(1,3))
p4000 <- ggplot(data = q4000, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p4000plus <- p4000 + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red")+
  labs(title="n=4000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

##

q8000 <- as.data.frame( rbind(cbind(Gn = 1/16, d = 3, q=qgum(n8000d3G16, 8000, 8000/16, 3)),
                              cbind(Gn = 1/16, d = 5, q=qgum(n8000d5G16, 8000, 8000/16, 5)),
                              cbind(Gn = 1/16, d = 7, q=qgum(n8000d7G16 , 8000, 8000/16, 7)),
                              cbind(Gn = 1/16, d = 10, q=qgum(n8000d10G16, 8000, 8000/16, 10)), 
                              cbind(Gn = 1/8, d = 3, q= qgum(n8000d3G8, 8000, 8000/8, 3 )),
                              cbind(Gn = 1/8, d = 5, q= qgum(n8000d5G8, 8000, 8000/8, 5 )),
                              cbind(Gn = 1/8, d = 7, q= qgum(n8000d7G8, 8000, 8000/8, 7)),
                              cbind(Gn = 1/8, d = 10, q= qgum(n8000d10G8, 8000, 8000/8, 10)),
                              cbind(Gn = 1/4, d = 3, q=qgum(n8000d3G4, 8000, 8000/4, 3)),
                              cbind(Gn = 1/4, d = 5, q=qgum(n8000d5G4, 8000, 8000/4, 5)),
                              cbind(Gn = 1/4, d = 7, q=qgum(n8000d7G4, 8000, 8000/4, 7)),
                              cbind(Gn = 1/4, d = 10, q=qgum(n8000d10G4, 8000, 8000/4, 10)) ) 
)
q8000$Gn <- factor(q8000$Gn)
q8000$d <- factor(q8000$d)

par(mfrow = c(1,3))
p8000 <- ggplot(data = q8000, aes(group=d, y=q)) + 
  geom_boxplot(aes(fill=d))
p8000plus <- p8000 + facet_wrap( ~ Gn) + geom_hline(yintercept=1, linetype="dashed", color = "red")+
  labs(title="n=8000", y= "Scaled Test Statistic" ) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

library(gridExtra)
grid.arrange(p1000plus, p2000plus, p4000plus, p8000plus)
