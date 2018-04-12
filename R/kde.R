
##' Suggest bandwith for local polynomial regression
##'
##' The bandwith is suggested coordinate wise to be
##' \deqn{0.9 \sqrt{5} \min\left(sd(x),\frac{IQR(x)}{1.349}\right) n ^{-\frac{1}{d+4}} (p+1)}
##' where p is the polynomial degree used and n is the number of coordinate points.
##' @param X A numeric matrix or vector of data coordinates
##' @param p Polynomial degree to base the suggestion on
##' @return a vector or scalar of suggested bandwiths
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats sd IQR
##' @export
suggestBandwith <- function(X, p){
    d <- ifelse(is.matrix(X),ncol(X),1)
    n <- ifelse(is.matrix(X),nrow(X),length(X))
    if(d == 1){
        bw <- 0.9 * sqrt(5) * min(stats::sd(X), stats::IQR(X) / 1.349) * n^(-(1/(d+4))) * (max(p,0) + 1)
    }else{
        bw <- 0.9 * sqrt(5) * apply(X,2,function(x){min(stats::sd(x),stats::IQR(x)/1.349)}) * n^(-(1/(d+4))) * (max(p,0) + 1)
    }
    return(bw)
}

##' Kernel Density Estimation
##'
##' Wrapper for the covafill reference class to do kernel density estimation.
##' @param X A numeric matrix or vector of data coordinates
##' @param bw Bandwith used
##' @param npred Number of coordinate wise equally spaced points at which the density is to be estimated. The numbers are repeated if the length is less than the dimension of the coordinates.
##' @param from Coordinate wise lower bound of points at which the density is to be estimated. The numbers are repeated if the length is less than the dimension of the coordinates.
##' @param to Coordinate wise upper bound of points at which the density is to be estimated. The numbers are repeated if the length is less than the dimension of the coordinates.
##' @return a list of coordinates and corresponding density estimates
##' @author Christoffer Moesgaard Albertsen
##' @export
kde <- function(X,bw = suggestBandwith(X,-1),npred = 100, from = min(X), to = max(X)){
   
    d <- ifelse(is.matrix(X),ncol(X),1)
    n <- ifelse(is.matrix(X),nrow(X),length(X))
    cf <- covafill(coord=X,obs=rep(1,n),p=-1L,h=bw)
    I <- d * (d + 2) * gamma(d/2) / (4 * pi ^ (d / 2))

    if(length(npred) < d)
        npred <- rep(npred, length = d)
    if(length(from) < d)
        from <- rep(from, length = d)
    if(length(to) < d)
        to <- rep(to, length = d)

    coords <- expand.grid(lapply(as.list(1:d),function(i)seq(from[i],to[i],length=npred[i])))
    dens <- apply(coords,1,function(x)I*cf$predict(matrix(x,1)))

    if(!is.null(colnames(X)))
        colnames(coords) <- colnames(X)
    
    return(list(coord=coords,density=dens))
}
