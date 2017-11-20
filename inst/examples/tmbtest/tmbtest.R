library(numDeriv)
library(TMB)
library(covafillr)

TMB::compile("tmbtest.cpp",CXXFLAGS=paste(cxxFlags(),"-std=c++14"))

dyn.load(dynlib("tmbtest"))




## Test one dim
coord <- as.matrix(expand.grid(seq(-1,1,0.001)))
ftrue <- function(x)sum(x^3)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

plot(coord,covObs)

dat <- list(coord = coord,
            covObs = covObs,
            p = 3,
            h = 0,
            d = 100)
dat$h <-  suggestBandwith(dat$coord,dat$p)

obj <- MakeADFun(data = dat,
                 parameters = list(x = c(0)),
                 DLL = "tmbtest")

fdiff <- Vectorize(function(x) obj$f(x) - x^3)
fdiff(seq(-1,1,0.01))


gdiff <- Vectorize(function(x) obj$gr(x) - 3*x^2)
gdiff(seq(-1,1,0.01))

cf <- covafill(dat$coord,dat$covObs,h=dat$h,p=as.integer(dat$p))
x0 <- seq(-1,1,0.1)
pr <- cf$predict(x0)

pr[,1]-x0^3
pr[,1]-2*x0^2


g1 <- 3*x0^2
g2 <- sapply(x0,obj$gr)
g3 <- pr[,2]

plot(x0,g1)
lines(x0,g2,col="red")
lines(x0,g3,col="blue")

## Test two dim

coord <- as.matrix(expand.grid(seq(-10,10,0.2),seq(-10,10,0.2)))
ftrue <- function(x)sum(x^3)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

dat <- list(coord = coord,
            covObs = covObs,
            p = 1,
            h = c(1,1))

obj <- MakeADFun(data = dat,
                 parameters = list(x = c(3.2,5.4)),
                 DLL = "tmbtest")


obj$fn(c(3.2,5.4))
obj$fn(c(0,0))
obj$fn(c(0,1))
obj$gr()
obj$he()


numDeriv::grad(obj$fn,c(0,1))
numDeriv::hessian(obj$fn,c(0,1))

fn <- Vectorize(function(x,y)obj$fn(c(x,y)))

system.time(fn(2,6))

x <- y <- seq(-5,5,0.1)

system.time(ztrue <- outer(x,y,Vectorize(function(x,y)ftrue(c(x,y)))))

system.time(zfit <- outer(x,y,fn))

lfit0 <- loess(covObs ~ x+y,data=data.frame(covObs=covObs,x=coord[,1],y=coord[,2]))
system.time(lfit <- outer(x,y,function(x,y)predict(lfit0,newdata=data.frame(x=x,y=y))))

par(mfrow=c(1,3))
image(x,y,ztrue)
image(x,y,zfit)
image(x,y,lfit)


rm(list=ls())
dyn.load(dynlib("tmbtest"))

coord <- as.matrix(expand.grid(seq(-10,10,0.5),
                               seq(-10,10,0.5),
                               seq(-10,10,0.5)))
ftrue <- function(x)sum(x^2)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

dat <- list(coord = coord,
            covObs = covObs,
            p = 2,
            h = c(1,1,1))

obj <- MakeADFun(data = dat,
                 parameters = list(x = c(3.2,5.4,0.2)),
                 DLL = "tmbtest")


###################
## One dimension ##
###################

dyn.load(dynlib("tmbtest"))

coord <- matrix(seq(-5,5,len=100),ncol=1)
ftrue <- function(x)x^2
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

library(ggplot2)

p <- 4

ggplot(data=data.frame(y=covObs,x=coord),aes(x=x,y=y)) +
    geom_point() +
    stat_covafill(polyDegree=p)

dat <- list(coord = coord,
            covObs = covObs,
            p = p,
            h = covafillr:::suggestBandwith(coord,p))

obj <- MakeADFun(data = dat,
                 parameters = list(x = c(2.5)),
                 DLL = "tmbtest")

obj$fn()
ftrue(2.5)
obj$gr()
2*2.5
obj$he(2.5)


ftrue <- function(x)x^3
coord <- matrix(seq(-5,5,len=150),ncol=1)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))
ggplot(data=data.frame(y=covObs,x=coord),aes(x=x,y=y)) +
    geom_point() +
    stat_covafill(polyDegree=p)

dat <- list(coord = coord,
            covObs = covObs,
            p = p,
            h = covafillr:::suggestBandwith(coord,p))

obj2 <- MakeADFun(data = dat,
                 parameters = list(x = c(2.5)),
                 DLL = "tmbtest")
obj2$fn()
ftrue(2.5)
obj2$gr()
3*2.5
obj2$he()

obj$he()
obj$fn()
obj$gr()


obj$env$Fun
obj2$env$Fun
