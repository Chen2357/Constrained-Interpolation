#' The Default Patching Dunction
#' 
#' The function used in the interpolations that use patching.
#' 
#' The percentage of \code{a} is given by \code{3*p^2-2*p^3}. Likewise, the percentage of \code{b} is given by \code{1-3*p^2+2*p^3}.
#' 
#' @param a A polynomial.
#' @param b A polynomial.
#' @param p A polynomial that describes the phase between \code{a} and \code{b}, with value between 0 and 1.
#' @return \code{(1-3*p^2+2*p^3)*(a-b)+b}
defaultPatchPolynomial <- function(a,b,p) {
    (((1 - (3 * (p ^ 2))) + (2 * (p ^ 3))) * (a - b)) + b
}

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

#' Quadratic Function Through Three Points
#' 
#' @param x A `pointData` type. It will be converted into `x` and `y`. Only the first three values will be used.
#' @param x The vector of x-coordinates. May be assigned directly.
#' @param y The vector of y-coordinates. May be assigned directly.
#' @return The quadratic polynomial that goes through \code{(x[1],y[1])}, \code{(x[2],y[2])}, and \code{(x[3],y[3])}.
quadraticPolynomial <- function(data, x = point.x(data), y = point.y(data)) {
    y[1] + ( ((y[2]-y[1])/(x[2]-x[1])) * polynomial(c(-x[1],10)) ) + ( (((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3])) * polynomial(c(-x[1],1)) * polynomial(c(-x[2],1)) )
}

#' Quadratic Function Through Two Points with Given Slope
#' 
#' This function produces the quadratic that goes through (x[1],y[1]) and (x[2],y[2]) and has a certain slope at (x[1],y[1])
#' 
#' @param x The vector of x-coordinates. Only the first two values will be used.
#' @param y The vector of y-coordinates. Only the first two values will be used.
#' @param s The slope at point \code{(x[1],y[1])}
#' @return The quadratic polynomial that goes through \code{(x[1],y[1])}, \code{(x[2],y[2])}, and has slope \code{s} at \code{(x[1],y[1])}.
projectQuadratic <- function(x,y,s) {
    y[1] + (s * polynomial(c(-x[1],1))) + ( ((y[2]-y[1]-s*(x[2]-x[1]))/((x[2]-x[1])^2)) * polynomial(c(x[1]^2,-2*x[1],1)) )
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

## ANCHOR Obsolete functions
#' Interpolation by Patching Quadratic Functions
#' 
#' Between any adjacent two points, the interpolated function is constructed by patching the two quadratic functions that go through the two points with the prescibed slope at each respective point.
#' 
#' If argument \code{slope} is missing, the interpolation will simply be patching two overlapping quadratic functions through two adjacent triplet of three points.
#' 
#' @param data A \code{pointData} type that stores all the points to be interpolated.
#' @param slope (Optional) The prescribed slopes at each points.
#' @param patch The function used for patching, uses \code{defaultPatchPolynomial} by default.
#' @return A piecewise polynomial.
interpolate.patchQuadratic <- function(data, slope, patch = defaultPatchPolynomial) {
    n <- length(data)
    leftBound <- point.x(data)[1:(n-1)]
    rightBound <- point.x(data)[2:n]
    polynomial <- vector(mode = "list", length = n-1)

    if (!missing(slope)) {
        for (i in 1:(n-1)) {
            polynomial[[i]] <- patch(
                projectQuadratic(point.x(data)[i:(i+1)],point.y(data)[i:(i+1)],slope[i]),
                projectQuadratic(point.x(data)[(i+1):i],point.y(data)[(i+1):i],slope[i+1]),
                percentagePolynomial(point.x(data)[i],point.x(data)[i+1])
            )
        }
    } else {
        # Original interpolation method
        polynomial[[1]] <- quadraticPolynomial(data[1:3])
        polynomial[[n-1]] <- quadraticPolynomial(data[n-2:n])
        for (i in 2:(n-2)) {
            polynomial[[i]] <- patch(
                quadraticPolynomial(data[(i-1):(i+1)]),
                quadraticPolynomial(data[i:(i+2)]),
                percentagePolynomial(point.x(data)[i],point.x(data)[i+1])
                )
        }
    }
    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

#' Interpolation by Patching Linear Functions
#' 
#' Between any adjacent two points, the interpolated function is constructed by patching the two linear functions that go through each respective point with the prescibed slope.
#' 
#' @param data A \code{pointData} type that stores all the points to be interpolated.
#' @param slope  The prescribed slopes at each points.
#' @param patch The function used for patching, uses \code{defaultPatchPolynomial} by default.
#' @return A piecewise polynomial.
interpolate.patchLinear <- function(data, slope, patch = defaultPatchPolynomial) {
    n <- length(data)
    leftBound <- point.x(data)[1:(n-1)]
    rightBound <- point.x(data)[2:n]
    polynomial <- vector(mode = "list", length = n-1)

    for (i in 1:(n-1)) {
        polynomial[[i]] <- patch(
            y[i] + (slope[i] * polynomial(-x[i],1)),
            y[i+1] + (slope[i+1] * polynomial(-x[i+1],1)),
            percentagePolynomial(point.x(data)[i],point.x(data)[i+1])
            )
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

#' Interpolation by Connecting Quadratic Functions
#' 
#' Between any adjacent two points, the interpolated function is constructed by connecting the two quadratic functions that (1) go through each respective point with the prescibed slope, (2) are smoothly connected, (3) have the same absolute curvature.
#' 
#' @param data A \code{pointData} type that stores all the points to be interpolated.
#' @param slope  The prescribed slopes at each points.
#' @return A piecewise polynomial.
interpolate.joinQuadratic <- function(data, slope) {
    n <- length(data)
    leftBound <- rep(NA, 2*n-2)
    rightBound <- rep(NA, 2*n-2)
    polynomial <- vector(mode = "list", length = 2*n-2)

    for (i in 1:(n-1)) {
        dx <- point.x(data)[i+1] - point.x(data)[i]
        a <- 2 * ((point.y(data)[i+1] - point.y(data)[i]) / dx - slope[i])
        ds <- slope[i+1] - slope[i]
        secondDerivative <- (a - ds + ifelse(a>ds,1,-1)*sqrt((a-ds)^2+ds^2)) / dx
        breakingPoint <- (point.x(data)[i+1] + point.x(data)[i] + ds / secondDerivative) / 2

        leftBound[2*i-1] <- point.x(data)[i]
        rightBound[2*i-1] <- breakingPoint
        polynomial[[2*i-1]] <- point.y(data)[i] + (slope[i] * polynomial(c(-point.x(data)[i],1))) + ((secondDerivative/2) * polynomial(c(point.x(data)[i]^2,-2*point.x(data)[i],1)))

        leftBound[2*i] <- breakingPoint
        rightBound[2*i] <- point.x(data)[i+1]
        polynomial[[2*i]] <- point.y(data)[i+1] + (slope[i+1] * polynomial(c(-point.x(data)[i+1],1))) - ((secondDerivative/2) * polynomial(c(point.x(data)[i+1]^2,-2*point.x(data)[i+1],1)))
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

## ANCHOR Flexible interpolation (New)
#' Interpolation by Patching Three Points Segments
#' 
#' @param data A `pointData` type that stores all the points to be interpolated.
#' @param solver Function that returns a polynomial that interpolate 3 points, must be in the form of `function(data)`.
#' @param patch The function used for patching, uses `defaultPatchPolynomial` by default. Must be in the form of `function(a,b,p)`.
#' @return A piecewise polynomial.
interpolate.patch.threePoint <- function(data, solver, patch = defaultPatchPolynomial) {
    if (length(data) < 3) stop("`data` must have length greater than 3")

    n <- length(data)
    leftBound <- point.x(data)[1:(n-1)]
    rightBound <- point.x(data)[2:n]
    polynomial <- vector(mode = "list", length = n-1)

    polynomial[[1]] <- solver(data[1:3])
    polynomial[[n-1]] <- solver(data[(n-2):n])
    previousPoly <- polynomial[[1]]
    for (i in 2:(n-2)) {
        thisPoly <- solver(data=data[i:(i+2)])
        polynomial[[i]] <- patch(previousPoly, thisPoly, percentagePolynomial(point.x(data)[i],point.x(data)[i+1]))
        previousPoly <- thisPoly
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

#' Interpolation by Patching Functions Generated at Each Point
#' 
#' @param data A `pointData` type that stores all the points to be interpolated.
#' @param slopes A `vector` type that stores all the slopes correspond to each point.
#' @param solver Function that returns a polynomial that is generated at a point with slope. Must be in the form of `function(data, slope)`.
#' @param patch The function used for patching, uses `defaultPatchPolynomial` by default. Must be in the form of `function(a,b,p)`.
#' @return A piecewise polynomial.
interpolate.patch.onePointSlope <- function(data, slopes, solver, patch = defaultPatchPolynomial) {
    n <- length(data)
    leftBound <- point.x(data)[1:(n-1)]
    rightBound <- point.x(data)[2:n]
    polynomial <- vector(mode = "list", length = n-1)

    previousPoly <- solver(data = data[1], slope = slopes[1])
    for (i in 1:(n-1)) {
        thisPoly <- solver(data = data[i+1], slope = slopes[i+1])
        polynomial[[i]] <- patch(previousPoly, thisPoly, percentagePolynomial(leftBound[i], rightBound[i]))
        previousPoly <- thisPoly
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

#' Interpolation by Joining Functions of Each segments
#' 
#' @param data A `pointData` type that stores all the points to be interpolated.
#' #' @param slopes A `vector` type that stores all the slopes correspond to each point.
#' @param solver Function that returns a polynomial that connect two points with respective slopes. Must be in the form of `function(data, slopes)`.
#' @param ... Additional parameters passed into `solver`.
#' @return A piecewise polynomial.
interpolate.twoPointSlope <- function(data, slopes, solver, ...) {
    n <- length(data)
    leftBound <- c()
    rightBound <- c()
    polynomial <- list()

    for (i in 1:(n-1)) {
        result <- solver(data = data[i:(i+1)], slopes = slopes[i:(i+1)], ...)
        if (class(result) == "piecewisePolynomial") {
            leftBound <- c(leftBound, result@leftBound)
            rightBound <- c(rightBound, result@rightBound)
            polynomial <- append(polynomial, result@polynomial)
        } else if (class(result) == "polynomial") {
            leftBound <- c(leftBound, point.x(data)[i])
            rightBound <- c(rightBound, point.x(data)[i+1])
            polynomial <- append(polynomial, result)
        } else if (class(result) == "numeric") {
            leftBound <- c(leftBound, point.x(data)[i])
            rightBound <- c(rightBound, point.x(data)[i+1])
            polynomial <- append(polynomial, polynomial(result))
        } else {
            stop(paste("Unrecognized return class of solver", class(result)))
        }
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}