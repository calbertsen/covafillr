
library(Matrix)
library(raster)

d <- 0.5
coords <- rbind(expand.grid(seq(1,5,d),seq(1,4,d)),
                expand.grid(seq(4,7,d),seq(5,8,d)))

tripletList <- array(dim=c(0,3))
delta <- 4
n <- dim(coords)[1]
for(i in 1:n)
    for(j in 1:n)
        if(sum((coords[i,]-coords[j,])^2)==d^2)
            tripletList <- rbind(tripletList,
                                 c(i,j,-1),
                                 c(i,i,1))
for(i in 1:n)
    tripletList <- rbind(tripletList,
                         c(i,i,delta))

Q <- sparseMatrix(i = tripletList[,1],
                  j = tripletList[,2],
                  x = tripletList[,3])

Sigma <- solve(Q)
CC <- cov2cor(Sigma)
image(cov2cor(Sigma))
CC[1:5,1:5]

Qsim <- Q + diag(-delta+0.001,dim(Q)[1],dim(Q)[2])
z <- mvtnorm::rmvnorm(1,sigma=as.matrix(solve(Qsim)))

r <- raster::raster(sp::SpatialPixelsDataFrame(coords,data=data.frame(z=as.vector(z))))
plot(r)

dat <- data.frame(z = as.vector(z),
                  x = coords[,1],
                  y = coords[,2])
pt0 <- unlist(coords[5,])

## fit <- lm(z ~ 1 + x + I(x^2) + y + I(y^2),
##           data=data.frame(z = dat$z,
##                           x = dat$x-pt0[1],
##                           y = dat$y-pt0[2]),
##           weights = Sigma[5,]/sum(Sigma[5,])
##           )
## z0 <- predict(fit,newdata=data.frame(x=0,y=0))

## B <- solve(t(chol(Q)))
## Ystar <- as.matrix(B)%*%dat$z

getVal <- Vectorize(function(x0,y0){
    indx <- which.min((dat$x-x0)^2 + (dat$y-y0)^2)
    if((dat$x[indx]-x0)^2 + (dat$y[indx]-y0)^2 > d^2)
        return(c(value = NA, grx = NA, gry=NA))
    
    X <- model.matrix(z ~ 1 + x + I(x^2) + y + I(y^2),
                      data=data.frame(z = dat$z,
                                      x = dat$x - x0,
                                      y = dat$y - y0))
    w <- CC[indx,] / sum(CC[indx,])
    Xstar <- w*X
    Ystar <- w*dat$z
    fit <- solve(t(Xstar)%*%Xstar)%*%t(Xstar) %*% Ystar
    c(value = fit[1], grx = fit[2], gry = fit[4])
})

newcoords <- expand.grid(seq(1,7,0.1),seq(1,8,0.1))
vals <- apply(newcoords,1,function(x)getVal(x[1],x[2]))

dev.new()
plot(r,col=rev(rainbow(99,start=0,end=1)),
     breaks=seq(min(minValue(r),minValue(r2)),max(maxValue(r),maxValue(r2)),length.out=100)
     ) 
dev.new()
r2 <- raster::raster(sp::SpatialPixelsDataFrame(newcoords,data=data.frame(z=as.vector(vals[1,]))))
plot(r2,col=rev(rainbow(99,start=0,end=1)),
     breaks=seq(min(minValue(r),minValue(r2)),max(maxValue(r),maxValue(r2)),length.out=100),
     main ="Smoothed field")
dev.new()
r3 <- raster::raster(sp::SpatialPixelsDataFrame(newcoords,data=data.frame(z=as.vector(vals[2,]))))
plot(r3, main = "x gradient of field")
dev.new()
r4 <- raster::raster(sp::SpatialPixelsDataFrame(newcoords,data=data.frame(z=as.vector(vals[3,]))))
plot(r4, main = "y gradient of field")



## Test with xtractomatic
## devtools::install_github("rmendels/xtractomatic")
library(xtractomatic)
library(raster)
library(rgdal)
library(Matrix)

d <- 0.5
llgrd <- expand.grid(seq(1,20,d),seq(53,59,d))+d/2
llgrd <- llgrd[abs(llgrd[,1]) < 179.995-d/2 & abs(llgrd[,2]) < 79.995-d/2,]
sst <- xtracto(llgrd[,1], llgrd[,2], rep('2015-11-21',dim(llgrd)[1]),"jplG1SST",xlen=0.001,ylen=0.001)
head(sst)

sst$mean[!is.finite(sst$mean)] <- NA
spdf <- SpatialPixelsDataFrame(as.matrix(llgrd),data=data.frame(z=unlist(sst$mean)))
rsst <- raster::raster(spdf)
proj4string(rsst) <- "+proj=longlat"
kom <- readOGR("/home/cmoe/kort/coastline","ne_10m_coastline")
kom <- spTransform(kom,CRS("+proj=longlat"))

plot(rsst)
plot(kom, main = "SST",add=TRUE, border ="black")


dat <- data.frame(z = sst$mean[is.finite(sst$mean)],
                  x = llgrd[is.finite(sst$mean),1],
                  y = llgrd[is.finite(sst$mean),2])
tripletList <- array(dim=c(0,3))
delta <- 0.005
n <- dim(dat)[1]
for(i in 1:n)
    for(j in 1:n)
        if(sum((dat[i,-1]-dat[j,-1])^2)==d^2)
            tripletList <- rbind(tripletList,
                                 c(i,j,-1),
                                 c(i,i,1))
for(i in 1:n)
    tripletList <- rbind(tripletList,
                         c(i,i,delta))
Q <- sparseMatrix(i = tripletList[,1],
                  j = tripletList[,2],
                  x = tripletList[,3])
Sigma <- solve(Q)
CC <- cov2cor(Sigma)

getVal <- Vectorize(function(x0,y0){
    indx <- which.min((dat$x-x0)^2 + (dat$y-y0)^2)
    if((dat$x[indx]-x0)^2 + (dat$y[indx]-y0)^2 > d^2)
        return(c(value = NA, grx = NA, gry=NA))
    X <- model.matrix(z ~ 1 + x + I(x^2) + y + I(y^2),
                      data=data.frame(z = dat$z,
                                      x = dat$x - x0,
                                      y = dat$y - y0))
    w <- CC[indx,] / sum(CC[indx,])
    Xstar <- w*X
    Ystar <- w*dat$z
    tryCatch({fit <- solve(t(Xstar)%*%Xstar)%*%t(Xstar) %*% Ystar
        return(c(value = fit[1], grx = fit[2], gry = fit[4]))
    },error=function(e)return(c(value = NA, grx = NA, gry=NA)))  
})

newcoords <- expand.grid(seq(min(llgrd[,1]),max(llgrd[,1]),0.1),
                         seq(min(llgrd[,2]),max(llgrd[,2]),0.1))
vals <- apply(newcoords,1,function(x)getVal(x[1],x[2]))

rssts <- raster::raster(SpatialPixelsDataFrame(as.matrix(newcoords),
                                               data=data.frame(z=unlist(vals[1,]))))

rp1 <- rsst
## rp1[rp1>20] <- NA
rp2 <- rssts
## rp2[rp2>20] <- NA
dev.new()
plot(rp1,col=rev(rainbow(99,start=0,end=1)),
     breaks=seq(min(minValue(rsst),minValue(rssts)),
                min(20,max(maxValue(rsst),maxValue(rssts))),length.out=100),
     main ="Observed")
plot(kom, main = "SST",add=TRUE, border ="black")
dev.new()
plot(rp2,col=rev(rainbow(99,start=0,end=1)),
     breaks=seq(min(minValue(rsst),minValue(rssts)),
                min(20,max(maxValue(rsst),maxValue(rssts))),length.out=100),
     main ="Smoothed field")
plot(kom, main = "SST",add=TRUE, border ="black")
