setGeneric("findSlope", function(finder,data) standardGeneric("findSlope"))

threePointQuadraticSlope <- function(data,z) {
    x = point.x(data)
    y = point.y(data)
    if (missing(z)) { z <- x[2] }

    return(
        (y[2]-y[1])/(x[2]-x[1]) + ((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3]) * (2*z - x[1] - x[2])
    )
}

quadraticSlopeFind <- function(data) {
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
    return(slopedPointData(x=point.x(data),y=point.y(data),slope=slope))
}

constrainedQuadratic.minSlope <- function(data) {
    x = point.x(data)
    y = point.y(data)

    # Solve for a(x-h)^2
    h <- (sqrt(y[1])*x[2] + sqrt(y[2]*x[1])) / (sqrt(y[1])+sqrt(y[2]))
    a <- (y[1])/((x[1]-h)^2)

    return(2*a*(x[1]-h))
}

constrainedQuadratic.maxSlope <- function(data) {
    x = point.x(data)
    y = point.y(data)

    h <- (sqrt(y[1])*x[2] + sqrt(y[2]*x[1])) / (sqrt(y[1])+sqrt(y[2]))
    if (is.infinite(h)) {
        return(0)
    }
    a <- (y[1])/((x[1]-h)^2)

    return(2*a*(x[2]-h))
}

constrainedQuadraticSlopeFind <- function(data) {
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
    return(slopedPointData(x=point.x(data),y=point.y(data),slope=slope))
}