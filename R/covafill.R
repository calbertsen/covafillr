#' A Reference Class for Local Polynomial Regression with covafill.
#' 
#' @field ptr External pointer to the covafill C++ object
#'
#' @examples
#' getRefClass('covafill')
#' fn <- function(x) x ^ 4 - x ^ 2
#' x <- runif(100,-3,3)
#' y <- fn(x) + rnorm(100,0,0.1)
#' cf <- covafill(coord = x,obs = y,h = 1.0,p = 5L)
#' cf$getDim()
#' cf$getDegree()
#' cf$getBandwith()
#' cf$setBandwith(1.0)
#' cf$getBandwith()
#' x0 <- seq(-1,1,0.1)
#' y0 <- cf$predict(x0)
#' par(mfrow=c(3,1))
#' plot(x0,y0[,1], main = "Function")
#' lines(x0,fn(x0))
#' plot(x0, y0[,2], main = "First derivative")
#' lines(x0, 4 * x0 ^ 3 - 2 * x0)
#' plot(x0, y0[,3], main = "Second derivative")
#' lines(x0, 3 * 4 * x0 ^ 2 - 2)
#' 
#' @export covafill
#' @importFrom methods setRefClass
#' @exportClass covafill
covafill <- setRefClass("covafill",
                        fields = list(ptr = "externalptr"),
                        methods = list(

                            initialize = function(.Object,
                                                  coord,
                                                  obs,
                                                  h = 1.0,
                                                  p = 2L,
                                                  ...){
                                "Method to initialize the covafill. coord is a matrix of coordinates, obs is a vector of corresponding observations, h is a vector of bandwiths, and p is the polynomial degree."
                                ## Check input
                                if(missing(coord) || missing(obs))
                                    stop("coord and obs must be specified")
                                if(!is.numvec(coord) & !is.nummat(coord))
                                    stop("coord must be a numeric vector or matrix.")

                                if(!is.numvec(obs))
                                    stop("Obs must be a numeric vector")
                                
                                if(!is.samelength(coord,obs))
                                    stop("obs and coord must have matching dimensions.")
                                if(is.numvec(coord))
                                    coord <- matrix(coord,ncol=1)
                                if(length(h) > dim(coord)[2])
                                    warning("Additional bandwiths are ignored.")

                                h <- rep(h,len = dim(coord)[2])

                                if(!is.numsca(p) & !is.numint(p))
                                    stop("p must be an integer.")

                                if(!is.numint(p))
                                    warning("p is coerced to an integer.")

                                p <- as.integer(p)

                                if(p < 1)
                                    stop("p must be 1 or greater")

                                ## Create pointer
                                ptr0 <- .Call("MakeFill",coord,obs,h,p,
                                              PACKAGE="covafillr")

                                initFields(ptr = ptr0)

                            },
                            predict = function(coord){
                                "Predict function value and derivatives with local polynomial regression at coord."
                                d <- .self$getDim()
                                if(!is.numvec(coord) & !is.nummat(coord))
                                    stop("coord must be a numeric vector or matrix.")
                                if(is.numvec(coord))
                                    coord <- matrix(coord,ncol=1)
                                if(dim(coord)[2] != d)
                                    stop(paste("coord must have",d,"columns."))

                                return(.Call("predictFill",.self$ptr,coord,
                                              PACKAGE="covafillr"))
                            },
                            residuals = function(excludeRadius){
                                "Get 'leave-neighborhood-out' residuals, i.e. local polynomial regression predictions excluding points within excludeRadius subtracted from the observation." 
                                if(!is.numsca(excludeRadius))
                                    stop("excludeRadius must be a numeric scalar")
                                return(.Call("lnoResiduals",.self$ptr,excludeRadius,
                                              PACKAGE="covafillr"))
                            },
                            getDim = function(){
                                "Get the dimension of the coordinates."
                                return(.Call("getFillDim",.self$ptr,
                                              PACKAGE="covafillr"))
                            },
                            getDegree = function(){
                                "Get the polynomial degree."
                                return(.Call("getFillDegree",.self$ptr,
                                              PACKAGE="covafillr"))
                            },
                            getBandwith = function(){
                                "Get the bandwith."
                                return(.Call("getFillBandwith",.self$ptr,
                                              PACKAGE="covafillr"))
                            },
                            setBandwith = function(h){
                                "Set the bandwith to h."
                                d <- .self$getDim()
                                if(length(h) > d)
                                    warning("Additional bandwiths are ignored.")
                                h <- rep(h,len = d)
                                return(.Call("setFillBandwith",.self$ptr,h,
                                              PACKAGE="covafillr"))
                            }
                            )
                        )
covafill$lock("ptr")
