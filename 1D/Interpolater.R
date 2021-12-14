#' Quadratic Function Through Two Points with Given Slope
#' 
#' This function produces the quadratic that goes through (x[1],y[1]) and (x[2],y[2]) and has a certain slope at (x[1],y[1])
#' 
#' @param x The vector of x-coordinates. Only the first two values will be used.
#' @param y The vector of y-coordinates. Only the first two values will be used.
#' @param s The slope at point \code{(x[1],y[1])}
#' @return The quadratic polynomial that goes through \code{(x[1],y[1])}, \code{(x[2],y[2])}, and has slope \code{s} at \code{(x[1],y[1])}.
projectQuadratic <- function(data, x = point.x(data), y = point.y(data), slope) {
    y[1] + (slope * polynomial(c(-x[1],1))) + ( ((y[2]-y[1]-slope*(x[2]-x[1]))/((x[2]-x[1])^2)) * polynomial(c(x[1]^2,-2*x[1],1)) )
}

patchProjectQuadratic <- function(data, slopes, patch = defaultPatchPolynomial) {
    patch(
        projectQuadratic(data[1:2], slope = slopes[1]),
        projectQuadratic(data[2:1], slope = slopes[2]),
        percentagePolynomial(point.x(data)[1], point.x(data)[2])
    )
}

#' Connecting  Two Quadratic Functions
#' 
#' Between two points, the function is constructed by connecting the two quadratic functions that (1) go through each respective point with the prescibed slope, (2) are smoothly connected, (3) have the same absolute curvature.
#' 
#' @param data A \code{pointData} type that stores all the points to be interpolated.
#' @param slope  The prescribed slopes at each points.
#' @return A piecewise polynomial.
joinQuadratic <- function(data, slope) {
    dx <- point.x(data)[2] - point.x(data)[1]
    a <- 2 * ((point.y(data)[2] - point.y(data)[1]) / dx - slope[1])
    ds <- slope[2] - slope[1]
    secondDerivative <- (a - ds + ifelse(a>ds,1,-1)*sqrt((a-ds)^2+ds^2)) / dx
    breakingPoint <- (point.x(data)[2] + point.x(data)[1] + ds / secondDerivative) / 2
    
    polynomial[[2*i-1]] <- point.y(data)[1] + (slope[1] * polynomial(c(-point.x(data)[1],1))) + ((secondDerivative/2) * polynomial(c(point.x(data)[1]^2,-2*point.x(data)[1],1)))
}

## ANCHOR Interpolation

#' Interpolation by Joining Functions of Each segments
#' 
#' @param data A `pointData` type that stores all the points to be interpolated.
#' @param slopes A `vector` type that stores all the slopes correspond to each point.
#' @param solver Function that returns a polynomial that connect two points with respective slopes. Must be in the form of `function(data, slopes)`.
#' @param ... Additional parameters passed into `solver`.
#' @return A piecewise polynomial.
interpolate.twoPointSlope <- function(data, slopes, solver, ...) {
    n <- length(data)
    x <- point.x(data)
    result <- piecewisePolynomial()

    for (i in 1:(n-1)) {
        poly <- solver(data = data[i:(i+1)], slopes = slopes[i:(i+1)], ...)
        result <- result %+% as.piecewisePolynomial(poly, x[i], x[i+1])
    }

    return(result)
}
