## Comparing sqrt(2logn) and D_n

nvals <- seq(1:10000)
Gvals <- nvals^.66
alpha <- 0.05
A <- 5
c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
a <- sqrt(2*log(nvals/Gvals)) #test transform multipliers
b5 <- 2*log(nvals/Gvals) + A/2 * log(log(nvals/Gvals)) - log(2/3 * gamma(A/2)) ##CORRECTED
D_5 <- (b5+c_alpha)/a #threshold
A <- 10
b10 <- 2*log(nvals/Gvals) + A/2 * log(log(nvals/Gvals)) - log(2/3 * gamma(A/2)) ##CORRECTED
D_10 <- (b10+c_alpha)/a #threshold
A <- 20
b20 <- 2*log(nvals/Gvals) + A/2 * log(log(nvals/Gvals)) - log(2/3 * gamma(A/2)) ##CORRECTED
D_20 <- (b20+c_alpha)/a #threshold
A <- 40
b40 <- 2*log(nvals/Gvals) + A/2 * log(log(nvals/Gvals)) - log(2/3 * gamma(A/2)) ##CORRECTED
D_40 <- (b40+c_alpha)/a #threshold

library(tidyverse)
df <- tibble(nvals, Limit = sqrt(2*log(nvals))+c_alpha/sqrt(2*log(nvals)),  D_5, D_10,  D_20,  D_40)
df1 <- df %>%
  gather(key = "variable", value = "value", -nvals)
ggplot(df1, aes(x = nvals, y = value)) + 
  geom_line(aes(color = variable)) + coord_cartesian(ylim = c(-15, 7)) + ggtitle("Transformed Critical Values", subtitle = "alpha = 0.05, G=n^(2/3)")

#ggplot(df, aes(x=nvals)) + geom_line(aes(y=Limit, color="black")) + geom_line(aes(y=D_40, color="darkblue")) + coord_cartesian(ylim = c(-30, 6))# + geom_line(aes(D_5, color="blue")) 
  