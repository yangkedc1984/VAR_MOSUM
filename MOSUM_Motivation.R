library(mosum)
x <- testData("custom", lengths = c(500,400,100,500), means = c(0,5,-5,0), sds = c(5,5,5,5))
ts.plot(x$x, ylab = "")
abline(v= c(500,900,1000), col = "red")

ms <- mosum(x$x,G=100, alpha = 0.01)
ms

plot(ms, display = "mosum")

library(wbs)
y <- testData("custom", lengths = c(100,100), means = c(0,5), sds = c(7,7))
bs <-wbs::sbs(y$x)
ts.plot(bs$x, ylab = "")
abline(v=100, col = "red")

res <- bs$res

ts.plot(abs(res[,4]), ylab = "")
abline(v=which.max(abs(res[,4])), col = "red")