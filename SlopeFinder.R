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

#' Slope Finding using Beta Method (3 points)
#' 
#' @param data A `pointData` type that stores the 3 points. It will be converted into `x` and `y`.
#' @param x A vector of x-coordinates of the 3 points. May be assigned directly.
#' @param y A vector of y-coordinates of the 3 points. May be assigned directly.
#' @param betaSolver A function to solve beta in the form of \code{function(A, B, b)} that returns a list containing beta. Can be `solve.beta` (default) or `solve.beta.cholesky`.
#' @return The slopes at the points of interest.
findSlope.beta.threePoints <- function(data, x = point.x(data), y = point.y(data), betaSolver = solve.beta) {
    delta21 <- x[2] - x[1]
    delta32 <- x[3] - x[2]

    # L_cluster matrix
    # d21^(-2)  0           -d21^(-2)   d21^(-1)    0           0
    # 0         0           d32^(-2)    0           -d32^(-2)   d32^(-1)
    # 0         0           0           0           1           0
    # 0         d21^(-1)    0           -d21^(-1)   0           0
    # 0         0           0           d32^(-1)    0           -d32^(-1)
    # 0         0           0           0           0           1
    L <- matrix(c(
        delta21^(-2),0,0,0,0,0,
        0,0,0,delta21^(-1),0,0,
        -delta21^(-2),delta32^(-2),0,0,0,0,
        delta21^(-1),0,0,-delta21^(-1),delta32^(-1),0,
        0,-delta32^(-2),1,0,0,0,
        0,delta32^(-1),0,0,-delta32^(-1),1),
        nrow = 6, ncol = 6
    )

    L_inv <- Solve(L)

    temp <- ifelse(y==0, 0, 1/(4*y))
    A <- t(L_inv) %*% diag(c(0, temp[1], 0, temp[2], 0, temp[3])) %*% L_inv
    
    temp <- ifelse(y==0, 1, 0)
    B <- diag(c(1, temp[1], 1, temp[2], 1, temp[3])) %*% L_inv
    b <- c(y[1], 0, y[2], 0, y[3], 0)

    beta <- betaSolver(A = A, B = B, b = b)
    p <- L_inv %*% beta

    return(p[c(2,4,6)])
}