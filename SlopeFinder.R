#' Slope of a Quadratic Through Three Points
#' 
#' @param data A \code{pointData} type that stores the 3 points that determine the quadratic function.
#' @param z The value at which the slope is evaluated.
#' @return The slope of the quadratic that goes through \code{data[1:3]} at value \code{z}.
threePointQuadraticSlope <- function(data,z) {
    x <- point.x(data)
    y <- point.y(data)
    if (missing(z)) z <- x[2]

    return(
        (y[2]-y[1])/(x[2]-x[1]) + ((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3]) * (2*z - x[1] - x[2])
    )
}

#' Slope Finding using Quadratics
#' 
#' @param data A \code{pointData} type.
#' @return For each point, the slope of the quadratic function that goes though it and the adjacent two points.
findSlope.quadratic <- function(data) {
    n <- length(data)
    slope <- rep(NA, n)
    for (i in 1:n) {
        if (i==1) {
            slope[1] <- threePointQuadraticSlope(data[1:3], point.x(data)[1])
        } else if (i==n) {
            slope[n] <- threePointQuadraticSlope(data[(n-2):n], point.x(data)[n])
        } else {
            slope[i] <- threePointQuadraticSlope(data[(i-1):(i+1)])
        }
    }
    return(slope)
}

#' Slope Lower Bound using Quadratic interpolation
#' 
#' @param data A \code{pointData} type that stores two points where the former point is the point of interest.
#' @return The slope at the point of interest such that the quadratic function going through the two points with the prescribed slope is nonnegative.
constrainedQuadratic.minSlope <- function(data) {
    x <- point.x(data)
    y <- point.y(data)

    if (y[1] == 0) {
        return(0)
    }

    # Solve for a(x-h)^2
    h <- (sqrt(y[1])*x[2] + sqrt(y[2]*x[1])) / (sqrt(y[1])+sqrt(y[2]))
    a <- (y[1])/((x[1]-h)^2)

    return(2*a*(x[1]-h))
}

#' Slope Upper Bound using Quadratic interpolation
#' 
#' @param data A \code{pointData} type that stores two points where the latter point is the point of interest.
#' @return The slope at the point of interest such that the quadratic function going through the two points with the prescribed slope is nonnegative.
constrainedQuadratic.maxSlope <- function(data) {
    x <- point.x(data)
    y <- point.y(data)

    if (y[2] == 0) {
        return(0)
    }

    h <- (sqrt(y[1])*x[2] + sqrt(y[2]*x[1])) / (sqrt(y[1])+sqrt(y[2]))
    if (is.infinite(h)) {
        return(0)
    }
    a <- (y[1])/((x[1]-h)^2)

    return(2*a*(x[2]-h))
}

#' Slope Finding using Quadratics with Constraint
#' 
#' @param data A \code{pointData} type.
#' @return For each point, the slope of the quadratic function that goes though it and the adjacent two points, but constrained such that the two quadratic functions extended from that point with the slope is always nonnegative.
findSlope.constrainedQuadratic <- function(data) {
    n <- length(data)
    slope <- rep(NA, n)
    for (i in 1:n) {
        if (i==1) {
            minSlope <- constrainedQuadratic.minSlope(data[1:2])
            slope[1] <- max(threePointQuadraticSlope(data[1:3], point.x(data)[1]), minSlope)
        } else if (i==n) {
            maxSlope <- constrainedQuadratic.minSlope(data[(n-1):n])
            slope[n] <- min(threePointQuadraticSlope(data[(n-2):n], point.x(data)[n]), maxSlope)
        } else {
            minSlope <- constrainedQuadratic.minSlope(data[i:(i+1)])
            maxSlope <- constrainedQuadratic.maxSlope(data[(i-1):i])
            slope[i] <- median(c(threePointQuadraticSlope(data[(i-1):(i+1)]), minSlope, maxSlope))
        }
    }
    return(slope)
}