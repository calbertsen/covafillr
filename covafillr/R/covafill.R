#' A Reference Class for Local Polynomial Regression with covafill.
#' 
#' @field ptr External pointer to the covafill C++ object
#'
#' @examples
#' getRefClass('covafill')
#' fn <- function(x) x ^ 4 - x ^ 2
#' x <- runif(2000,-3,3)
#' y <- fn(x) + rnorm(2000,0,0.1)
#' cf <- covafill(coord = x,obs = y,p = 5L)
#' cf$getDim()
#' cf$getDegree()
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
#' cf$setBandwith(1.0)
#' cf$getBandwith()
#' 
#' @export covafill
#' @importFrom methods setRefClass new 
#' @exportClass covafill
#' @useDynLib covafillr, .registration = TRUE, .fixes = "C_"
covafill <- setRefClass("covafill",
                        fields = list(ptr = "externalptr"),
                        methods = list(

                            initialize = function(coord,
                                                  obs,
                                                  h = suggestBandwith(coord,p),
                                                  p = 3L,
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
                                if(is.numvec(coord)){
                                    coord <- matrix(coord,ncol=1)
                                }
                                if(length(h) > dim(coord)[2])
                                    warning("Additional bandwiths are ignored.")

                                h <- rep(h,len = dim(coord)[2])

                                if(!is.numsca(p) & !is.numint(p))
                                    stop("p must be an integer.")

                                if(!is.numint(p))
                                    warning("p is coerced to an integer.")

                                p <- as.integer(p)

                                if(p < -1)
                                    stop("p must be -1 or greater")

                                ## Create pointer
                                ptr0 <- .Call(C_MakeFill,coord,obs,h,p)

                                initFields(ptr = ptr0)

                            },
                            copy = function(shallow=FALSE){
                                ##stop("A covafill object can not be copied.")
                                if(shallow)
                                    stop("covafill only allows non-shallow copy")
                                obs <- .Call(C_getFillObservations,.self$ptr)
                                coord <- .Call(C_getFillCoordinates,.self$ptr)
                                bw <- .self$getBandwith()
                                p <- .self$getDegree()
                                do.call("covafill",list(h=bw,p=p,obs=obs,coord=coord))
                            },
                            predict = function(coord, se.fit=FALSE){
                                "Predict function value and derivatives with local polynomial regression at coord. If se.fit=TRUE a list is returned with estimates and their standard deviations."
                                d <- .self$getDim()
                                p <- .self$getDegree()
                                
                                if(!is.numvec(coord) & !is.nummat(coord))
                                    stop("coord must be a numeric vector or matrix.")
                                if(is.numvec(coord))
                                    coord <- matrix(coord,ncol=1)
                                if(dim(coord)[2] != d)
                                    stop(paste("coord must have",d,"columns."))

                                if(se.fit){
                                    val <- .Call(C_predictFillSE,.self$ptr,coord)
                                    val[[2]] <- sqrt(val[[2]])
                                }else{
                                    val <- .Call(C_predictFill,.self$ptr,coord)
                                }

                                nest <- ifelse(se.fit,dim(val[[1]])[2],dim(val)[2])
                                
                                if(is.null(colnames(coord))){
                                    cnam <- 1:d
                                }else{
                                    cnam <- colnames(coord)
                                }

                                cnamfin <- character(nest)
                                cnamfin[1] <- 'fn'
                                
                                if(p >= 1){
                                    cnamfin[2:(1+d)] <- paste('gr',cnam,sep='_')

                                    if(p >= 2){
                                        cn2 <- unlist(sapply(1:d,function(i)
                                            paste(cnam[i],cnam[i:d],sep='_')))
                                        cnamfin[(2+d):(2+d+length(cn2)-1)] <- paste('gr',
                                                                                    cn2,
                                                                                    sep='_')
                                    }
                                    if(p > 2){
                                        cnK <- as.vector(sapply(3:p,
                                                                function(k)
                                                                    unlist(sapply(1:d,
                                                                                  function(j)
                                                                                      paste(rep(cnam[j],k),collapse='_')
                                                                                  ))
                                                                ))
                                        cnamfin[tail(1:length(cnamfin),
                                                     length(cnK))] <- paste('gr',cnK,sep='_')
                                    }
                                }

                                if(is.null(rownames(coord))){
                                    rnam <- 1:dim(coord)[1]
                                }else{
                                    rnam <- rownames(coord)
                                }

                                if(se.fit){
                                    names(val) <- c("fit","se.fit")
                                    colnames(val[[1]]) <- colnames(val[[2]]) <- cnamfin
                                    rownames(val[[1]]) <- rownames(val[[2]]) <- rnam

                                }else{
                                    colnames(val) <- cnamfin
                                    rownames(val) <- rnam
                                }
                                
                                return(val)
                            },
                            residuals = function(excludeRadius){
                                "Get 'leave-neighborhood-out' residuals, i.e. local polynomial regression predictions excluding points within excludeRadius subtracted from the observation." 
                                if(!is.numsca(excludeRadius))
                                    stop("excludeRadius must be a numeric scalar")
                                return(.Call(C_lnoResiduals,.self$ptr,excludeRadius))
                            },
                            getDim = function(){
                                "Get the dimension of the coordinates."
                                return(.Call(C_getFillDim,.self$ptr))
                            },
                            getDegree = function(){
                                "Get the polynomial degree."
                                return(.Call(C_getFillDegree,.self$ptr))
                            },
                            getBandwith = function(){
                                "Get the bandwith."
                                return(.Call(C_getFillBandwith,.self$ptr))
                            },
                            setBandwith = function(h){
                                "Set the bandwith to h."
                                d <- .self$getDim()
                                if(length(h) > d)
                                    warning("Additional bandwiths are ignored.")
                                h <- rep(h,len = d)
                                return(.Call(C_setFillBandwith,.self$ptr,h))
                            },
                            show = function(){
                                cat("\ncovafill Local Polynomial Regression Object\n\n")
                                obs <- .Call(C_getFillObservations,.self$ptr)
                                coord <- .Call(C_getFillCoordinates,.self$ptr)
                                bw <- .self$getBandwith()
                                p <- .self$getDegree()
                                cat("\nPolynomial degree:",p,"\n")
                                cat("\nBandwith:",bw,"\n")
                                cat("\nCoordinates used:\n")
                                print(summary(coord))
                                cat("\nObservations used:\n")
                                print(summary(obs))
                                invisible(list(h=bw,p=p,obs=obs,coord=coord))
                            }
                        )
                        )
covafill$lock("ptr")
