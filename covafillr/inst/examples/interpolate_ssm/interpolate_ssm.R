library(TMB)
library(covafillr)

compile("interpolate_ssm.cpp",CXXFLAGS = paste(cxxFlags(),"-std=c++14"))

dyn.load(dynlib("interpolate_ssm"))


coord <- as.matrix(expand.grid(seq(-1,1,len=5),seq(-1,1,len=5)))
ftrue <- function(x) sum(diag(1,2) %*% (0.5 * x * x - c(1,-1) * x))
covObsTrue <- apply(coord,1,function(x)ftrue(x) )##+ rnorm(1,0,0.0))

n <- 1000
dt <- 1/100
Xt <- matrix(0,2,n*(1/dt))
Xt[,1] <- 0.5
for(i in 2:ncol(Xt)){
    v <- numDeriv::grad(ftrue,Xt[,i-1])
    Xt[,i] <- Xt[,i-1] + -v*dt + rnorm(2,0,0.1 * sqrt(dt))
}
X <- Xt[,(1:ncol(Xt))*dt == round((1:ncol(Xt))*dt)]
plot(t(X),type="l")
plot(X[1,],type="l")

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
