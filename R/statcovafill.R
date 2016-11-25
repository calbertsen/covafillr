
#' @importFrom ggplot2 ggproto Stat layer ggproto
#'
#' @keywords internal
StatCovafill <- ggplot2::ggproto("StatCovafill", ggplot2::Stat,
                        required_as = c("x", "y"),
                        compute_group = function(data, scales, params,
                                                 n,
                                                 bandwith,
                                                 polyDegree,
                                                 level,
                                                 se){

                            if(is.null(bandwith))
                                bandwith <- suggestBandwith(data$x,polyDegree)

                            cf <- covafill(coord = data$x,
                                           obs = data$y,
                                           h = bandwith,
                                           p = as.integer(polyDegree))

                            rng <- range(data$x, na.rm = TRUE)
                            xseq <-seq(rng[1],rng[2],length=n)
                            pred <- cf$predict(xseq,se.fit = se)

                            if(se){
                                fac <- qnorm(level)
                                return(data.frame(x = xseq,
                                                  y = unname(pred$fit[,1]),
                                                  ymin = unname(pred$fit[,1] - fac * pred$se.fit[,1]),
                                                  ymax = unname(pred$fit[,1] + fac * pred$se.fit[,1]),
                                                  se = unname(pred$se.fit[,1])))
                            }else{
                                return(data.frame(x = xseq,
                                                  y = unname(pred[,1])))
                            }
                        }
                        )


##' Add a covafill smoother to an (x,y) plot
##'
##' As an extention to the \code{ggplot2} package, the function adds a covafill fit to an (x,y) plot. The fit is predicted to points on the interval range(x).
##'  
##' @param mapping Set of mappings created by 'aes' from the \code{ggplot2} package. The same as \code{ggplot2::stat_smooth}.
##' @param data The data to be displayed in this layer. The same as \code{ggplot2::stat_smooth}.
##' @param geom The same as \code{ggplot2::stat_smooth}.
##' @param position Position adjustments. The same as \code{ggplot2::stat_smooth}.
##' @param na.rm Not used
##' @param show.legend Should this legend be displayed? The same as \code{ggplot2::stat_smooth}.
##' @param inherit.aes The same as \code{ggplot2::stat_smooth}.
##' @param n Number of points to do prediction on.
##' @param bandwith Bandwith used in covafill. Uses \code{suggestBandwith} by default.
##' @param polyDegree Polynomial degree to use in covafill.
##' @param level Level of confidence interval to use.
##' @param se Should confidence intervals be displayed?
##' @param ... Other arguments passed to \code{layer}.
##' @return A \code{ggplot2} \code{layer}.
##' @author Christoffer Moesgaard Albertsen
##'
##' @seealso \code{\link[ggplot2]{stat_smooth}}
##' 
##' @export
stat_covafill <- function(mapping = NULL, data = NULL, geom = "smooth",
                    position = "identity", na.rm = FALSE, show.legend = NA, 
                    inherit.aes = TRUE, n = 50, bandwith = NULL, polyDegree = 3L, level = 0.95, se = TRUE, 
                    ...) {
    requireNamespace("ggplot2",quietly = TRUE)
    ggplot2::layer(
        stat = StatCovafill, data = data, mapping = mapping, geom = geom, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(n = n, bandwith = bandwith, polyDegree = polyDegree, level=level, se = se, ...)
    )
}


