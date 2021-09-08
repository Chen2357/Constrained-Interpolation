#' Percentage Polynomial
#' 
#' The linear function that is 0 when evaluated at \code{min} and 1 when evaluated at \code{max}
#' 
#' @param min The value at which the function is 0.
#' @param max The value at which the function is 1.
#' @return \code{(x-min)/(max-min)}
percentagePolynomial <- function(min,max) {
    (1/(max-min)) * polynomial(c(-min,1))
}

#' Linear Polynomial Given Point and Slope
#' 
#' @param x A `pointData` type. It will be converted into `x` and `y`. Only the first three values will be used.
#' @param x The vector of x-coordinates. May be assigned directly.
#' @param y The vector of y-coordinates. May be assigned directly.
#' @param slope The slope at the point.
#' @return The linear polynomial that goes through the point with the given slope.
linearPolynomial <- function(data, x = point.x(data), y = point.y(data), slope) {
    if (missing(slope)) {
        slope <- (y[2] - y[1]) / (x[2] - x[1])
    }
    return(polynomial(c(y[1] - slope * x[1], slope)))
}

#' Quadratic Function Through Three Points
#' 
#' @param x A `pointData` type. It will be converted into `x` and `y`. Only the first three values will be used.
#' @param x The vector of x-coordinates. May be assigned directly.
#' @param y The vector of y-coordinates. May be assigned directly.
#' @return The quadratic polynomial that goes through \code{(x[1],y[1])}, \code{(x[2],y[2])}, and \code{(x[3],y[3])}.
quadraticPolynomial <- function(data, x = point.x(data), y = point.y(data)) {
    y[1] + ( ((y[2]-y[1])/(x[2]-x[1])) * polynomial(c(-x[1],1)) ) + ( (((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3])) * polynomial(c(-x[1],1)) * polynomial(c(-x[2],1)) )
}

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

#' Quadratic Function Through a Point and with a Given Slope and with Given y-value of Extrema
#' 
#' @param data A `pointData` type containing a signle point. It will be converted to `x` and `y`
#' @param x The x-coordinate of a point. May be assigned directly.
#' @param y The y-coordinate of a point. May be assigned directly.
#' @param slope The slope at point `(x,y)``
#' @param extreme The y-value of the extrema.
#' @return The quadratic that goes through `(x,y)` with slope `k` and tangent to the `y=extrema`.
quadratic.point.slope.extrema <- function(data, x = point.x(data), y = point.y(data), slope, extrema = 0, tol = sqrt(.Machine$double.eps)) {
    if(abs(y-extrema) < tol) return(polynomial(y))
    slope^2/(4*(y-extrema))*polynomial(c(x^2,-2*x,1)) + slope * polynomial(c(-x,1)) + y
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

restrictedRange <- function(data, slope, tau, patch = patch.fifthDegree) {
    x0 <- point.x(data)
    k <- slope
    b <- point.y(data)

    p <- polynomial(c(b-k*x0, k))

    mu <- k^2 / (tau - abs(b))
    delta <- (tau - abs(b)) / k

    if (delta >= 1) {
        result <- patching(
            list(
                polynomial(0),
                p,
                polynomial(0)
            ),
            x0 + c(-1,0,1),
            patch
        )
    } else {
        q <- polynomial(c(0,0,mu/4))
        result <- patching(
            list(
                polynomial(-tau),
                p+q,
                p,
                p-q,
                polynomial(tau)
            ), 
            x0 + delta * c(-2*sqrt(2),-1,0, 1,2*sqrt(2)),
            patch
        )
    }
    return(result)
}

## ANCHOR Interpolation

#' Interpolation by Patching Three Points Segments
#' 
#' @param data A `pointData` type that stores all the points to be interpolated.
#' @param solver Function that returns a polynomial that interpolate 3 points, must be in the form of `function(data)`.
#' @param patch The function used for patching, uses `defaultPatchPolynomial` by default. Must be in the form of `function(a,b,p)`.
#' @return A piecewise polynomial.
interpolate.patch.threePoint <- function(data, solver, patch = patch.fifthDegree) {
    if (length(data) < 3) stop("`data` must have length of at least 3")

    n <- length(data)
    x <- point.x(data)

    func <- list()
    breaks <- x[2:(n-1)]

    for (i in 2:(n-1)) {
        func <- c(func, solver(data=data[(i-1):(i+1)]))
    }

    return(patching(func, breaks, patch))
}

#' Interpolation by Patching Functions Generated at Each Point
#' 
#' @param data A `pointData` type that stores all the points to be interpolated.
#' @param slopes A `vector` type that stores all the slopes correspond to each point.
#' @param solver Function that returns a polynomial that is generated at a point with slope. Must be in the form of `function(data, slope)`.
#' @param patch The function used for patching, uses `defaultPatchPolynomial` by default. Must be in the form of `function(a,b,p)`.
#' @return A piecewise polynomial.
interpolate.patch.onePointSlope <- function(data, slopes, solver, patch = patch.fifthDegree) {
    n <- length(data)
    x <- point.x(data)
    result <- piecewisePolynomial()

    func <- list()

    for (i in seq_len(n)) {
        func <- c(func, solver(data = data[i], slope = slopes[i]))
    }

    return(patching(func, x, patch))
}

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