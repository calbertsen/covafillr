library(numDeriv)
library(TMB)
library(covafillr)

system("touch tmbtest.cpp")
TMB::compile("tmbtest.cpp",CXXFLAGS=paste(cxxFlags(),"-std=c++14"))

dyn.load(dynlib("tmbtest"))

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
