library(ggplot2)
library(covafillr)

x <- seq(-2,2,len=300)
pol <- function(x)1 + 0.5 * x - 1.5 * x^2 + 2 * x ^ 3
y <- rnorm(length(x),pol(x), 2)

data <- data.frame(x=x,y=y)

cf <- covafill(coord = data$x,
               obs = data$y,
               h = 2.0,
               p = as.integer(3))


cf$predict(0,TRUE)

ggplot(aes(x=x,y=y),data=data) +
    geom_point() +
    stat_covafill(bandwith=1,polyDegree=3L) +
    geom_line(aes(x=x,y=y),data=data.frame(x=seq(-2,2,len=1000),y=pol(seq(-2,2,len=1000))),col="red")
