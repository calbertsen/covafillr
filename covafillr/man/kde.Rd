% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kde.R
\name{kde}
\alias{kde}
\title{Kernel Density Estimation}
\usage{
kde(X, bw = suggestBandwith(X, -1), npred = 100, from = min(X), to = max(X))
}
\arguments{
\item{X}{A numeric matrix or vector of data coordinates}

\item{bw}{Bandwith used}

\item{npred}{Number of coordinate wise equally spaced points at which the density is to be estimated. The numbers are repeated if the length is less than the dimension of the coordinates.}

\item{from}{Coordinate wise lower bound of points at which the density is to be estimated. The numbers are repeated if the length is less than the dimension of the coordinates.}

\item{to}{Coordinate wise upper bound of points at which the density is to be estimated. The numbers are repeated if the length is less than the dimension of the coordinates.}
}
\value{
a list of coordinates and corresponding density estimates
}
\description{
Wrapper for the covafill reference class to do kernel density estimation.
}
\author{
Christoffer Moesgaard Albertsen
}
