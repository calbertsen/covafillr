library(covafillr)

coord <- as.matrix(expand.grid(seq(-10,10,0.2),seq(-10,10,0.2)))
ftrue <- function(x)sum(x^3) + prod(x)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

## Test smoothField
cf <- covafill(coord=coord,obs=covObs,h=c(1,1),p=2L)
cf2 <- covafill(coord=coord,obs=covObs)

cf$accessors()



val1 <- cf$predict(matrix(c(0,0),1,2))
xx <- matrix(runif(200,-10,10),100,2)
val2 <- cf$predict(xx)

## Test smoothTree

ct <- covatree(coord=coord,obs=covObs,h=c(1,1),p=2L, minLeft = 100)
val3 <- ct$predict(xx)


## Cross validate - this takes a lot of time!

bw <- seq(0.5,5,0.5)
cv <- numeric(length(bw))
for(i in 1:length(bw)){
    cf$setBandwith(bw[i])
    res <- cf$residuals(0.1)
    cv[i] <- mean(res ^ 2)
}
    

plot(bw,cv)
