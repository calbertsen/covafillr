library(TMB)
library(covafillr)

compile("interpolate.cpp",CXXFLAGS = paste(cxxFlags(),"-std=c++14"))

dyn.load(dynlib("interpolate"))

coord <- as.matrix(expand.grid(seq(-1,1,len=5),seq(-1,1,len=5)))
ftrue <- function(x)sum(x^2)
covObsTrue <- apply(coord,1,function(x)ftrue(x) )##+ rnorm(1,0,0.0))

n <- 5000
X <- rbind(runif(n,-1,1),
           runif(n,-1,1))
obs <- apply(X,2,ftrue) + rnorm(n,0,0.1)

dat <- list(coord = coord,
            p = 2,
            h = 0,
            d = 2,
            obs = obs,
            x = X,
            repx = seq(-1,1,0.05),
            repy = seq(-1,1,0.05))
dat$h <- covafillr:::suggestBandwith(coord,dat$p)
pars <- list(covObs = rep(0,nrow(coord)),
             logSd = 0)

obj <- MakeADFun(data = dat,
                 parameters = pars,
                 DLL = "interpolate")

opt <- nlminb(obj$par,obj$fn,obj$gr,obj$he,control=list(iter.max=1000,eval.max=1000))

opt


rp <- obj$report()

contour(dat$repx,dat$repy,rp$predXY)
contour(dat$repx,dat$repy,outer(dat$repx,dat$repy,Vectorize(function(x,y)ftrue(c(x,y)))),col="red",add=TRUE)
