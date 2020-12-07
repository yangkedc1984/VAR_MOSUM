library(mosum)
library(tidyverse)
## Detecting the COVID-19 dip

ftse_28 <- readxl::read_xlsx("ftse100_28_02_2020.xlsx", sheet = 3)
ftse28.new <-  tibble(FTSE100 = rev(ftse_28$`FTSE 100`))
ftse_28$Date <- as.Date( ftse_28$Date, '%Y-%m-%d')
ftse28.new$Date <- rev(ftse_28$Date)
head(ftse28.new)

ggplot( data = ftse28.new, aes(Date, FTSE100 )) + geom_line() 

ms <- mosum(x = ftse28.new$FTSE100, G = 10, alpha = 0.1 ) #calc mosum
ms


plot(ms, xlab = "Date")
title(main = "FTSE100 in the time of COVID-19")
axis(side = 1, at = ftse28.new$Date, labels = ftse28.new$Date)

#restrict to - weeks
msd <- mosum(x = ftse28.new$FTSE100[70:129], G = 10, alpha = 0.1 ) #calc mosum
msd
plot(msd, xlab = "Date"  )

#ftse28.new$mosum <- ms$rollsums
#ggplot() + geom_line(data = ftse28.new, aes(Date, `FTSE 100` ), color = "black") + geom_line(data = ftse28.new, aes(x= Date, y = mosum), color = "red") +
#   geom_vline(xintercept = ms$cpts[1], colour = "blue")


## Diff data
ftse_diff <- diff(ftse28.new$`FTSE 100`) #diff
df <- tibble(Date = ftse28.new$Date[1:129], Diff = ftse_diff) #tibble
head(df)
ggplot( data = df, aes(Date, Diff )) + geom_line() 

msd <- mosum(x = df$Diff, G = 15, alpha = 0.1 )
#############################################################################
#USD - CNY
library(lubridate)
exch <- readxl::read_xlsx("USD_CNY.xlsx")
#exch$Date <- as.Date( exch$Date, '%m_%d,%Y')
exch$Date <- mdy(exch$Date)
exch_price <- exch[,1:2]
exch_price <- exch_price[order(exch_price$Date),] #reverse
head(exch_price)

ggplot( data = exch_price, aes(Date, Price )) + geom_line() 

ms <- mosum(x = exch_price$Price, G = 20, alpha = 0.1 ) #calc mosum
ms

plot(ms, xlab = "Date")
#title(main = "USD - CNY, 2019/03/18 - 2020/03/18")
#axis(side = 1, at = ftse28.new$Date, labels = ftse28.new$Date)
