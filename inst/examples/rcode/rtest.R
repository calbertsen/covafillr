library(covafillr)


## Dim: 1

X <- rnorm(1000,0,1)
Y <- cos(X)+rnorm(1000,0,0.1)

plot(X,Y)

plot(density(X,kernel="epanechnikov"))
cf <- covafill(coord=X,obs=rep(1,length(X)),p=-1L)
cf$predict(0,TRUE)
#I <- sum(cf$predict(seq(min(X),max(X),0.001)))*0.001
I <- 1/(1*(1+2)*gamma(0.5)/(4*pi^0.5))
lines(seq(-3,3,0.1),cf$predict(seq(-3,3,0.1))[,1]/I,col="red")

test <- kde(X)
lines(test[[1]][,1],test[[2]],col="blue")

cf$predict(0,TRUE)
cf$predict(0,FALSE)



library(ggplot2)

ggplot(data=data.frame(X=X,Y=Y),aes(x=X,y=Y)) +
    geom_point() +
    stat_covafill(polyDegree=6) +
    geom_line(aes(x=X,y=Y),data=data.frame(X=seq(min(X),max(X),len=1000),Y=cos(seq(min(X),max(X),len=1000))),col="red")

#Dim: 2

set.seed(123)
X <- t(replicate(1000,as.vector(t(chol(matrix(c(1,0.9,0.9,1),2))) %*% rnorm(2))))
Y <- rowSums(cos(X))


zkde <- covafillr:::kde(X,npred=25)

x <- sort(unique(zkde$coord[,1]))
y <- sort(unique(zkde$coord[,2]))
z <- matrix(zkde$density,length(x),length(y))

contour(x,y,z)

dev.new()
contour((z2 <- MASS::kde2d(X[,1],X[,2])))

dev.new()
contour(z2$z-z)

x <- y <- seq(min(X),max(X),len=100)


cf <- covafill(coord=X,obs=Y,p=-1L,h=c(1,1))
ff <- Vectorize(function(x,y)cf$predict(cbind(x,y))[1])
z <- outer(x,y,ff)

contour(x,y,z)


cf$predict(cbind(0,0))

cf <- covafill(X,Y,p=3L)
ct <- covatree(X,Y,p=5L)
cf$predict(matrix(0,1,2))    
ct$predict(matrix(0.1,1,2))



