defaultPatchPolynomial <- function(a,b,p) {
    (((1 - (3 * (p ^ 2))) + (2 * (p ^ 3))) * (a - b)) + b
}

percentagePolynomial <- function(min,max) {
    (1/(max-min)) * polynomial(c(-min,1))
}

quadraticPolynomial <- function(x,y) {
    y[1] + ( ((y[2]-y[1])/(x[2]-x[1])) * polynomial(c(-x[1],10)) ) + ( (((y[2]-y[1])/(x[2]-x[1])-(y[3]-y[2])/(x[3]-x[2]))/(x[1]-x[3])) * polynomial(c(-x[1],1)) * polynomial(c(-x[2],1)) )
}

# produces the quadratic that goes through (x[1],y[1]) and (x[2],y[2]) and has a certain slope at (x[1],y[1])
projectQuadratic <- function(x,y,s) {
    y[1] + (s * polynomial(c(-x[1],1))) + ( ((y[2]-y[1]-s*(x[2]-x[1]))/((x[2]-x[1])^2)) * polynomial(c(x[1]^2,-2*x[1],1)) )
}

quadraticPatchInterpolate <- function(data, patch = defaultPatchPolynomial) {
    n <- length(data)
    leftBound <- point.x(data)[1:(n-1)]
    rightBound <- point.x(data)[2:n]
    polynomial <- vector(mode = "list", length = n-1)

    if (is(data,"slopedPointData")) {
        for (i in 1:(n-1)) {
            polynomial[[i]] = patch(
                projectQuadratic(point.x(data)[i:(i+1)],point.y(data)[i:(i+1)],slope(data)[i]),
                projectQuadratic(point.x(data)[(i+1):i],point.y(data)[(i+1):i],slope(data)[i+1]),
                percentagePolynomial(point.x(data)[i],point.x(data)[i+1])
                )
            # validObject(polynomial[[i]])
        }
    } else if (is(data,"pointData")) {
        # Original interpolation method
        polynomial[[1]] = quadraticPolynomial(point.x(data)[1:3],y(data)[1:3])
        polynomial[[n-1]] = quadraticPolynomial(point.x(data)[n-2:n],y(data)[n-2:n])
        for (i in 2:(n-2)) {
            polynomial[[i]] = patch(
                quadraticPolynomial(point.x(data)[(i-1):(i+1)],point.y(data)[(i-1):(i+1)]),
                quadraticPolynomial(point.x(data)[i:(i+2)],point.y(data)[i:(i+2)]),
                percentagePolynomial(point.x(data)[i],point.x(data)[i+1])
                )
        }
    } else {
        warning("Unsupported data class")
        return(NULL)
    }
    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

linearPatchInterpolate <- function(data, patch = defaultPatchPolynomial) {
    n <- length(data)
    leftBound <- point.x(data)[1:(n-1)]
    rightBound <- point.x(data)[2:n]
    polynomial <- vector(mode = "list", length = n-1)

    for (i in 1:(n-1)) {
        polynomial[[i]] = patch(
            y[i] + (slope(data)[i] * polynomial(-x[i],1)),
            y[i+1] + (slope(data)[i+1] * polynomial(-x[i+1],1)),
            percentagePolynomial(point.x(data)[i],point.x(data)[i+1])
            )
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

quadraticJointInterpolate <- function(data) {
    n <- length(data)
    leftBound <- rep(NA, 2*n-2)
    rightBound <- rep(NA, 2*n-2)
    polynomial <- vector(mode = "list", length = 2*n-2)

    for (i in 1:(n-1)) {
        dx <- point.x(data)[i+1] - point.x(data)[i]
        a <- 2 * ((point.y(data)[i+1] - point.y(data)[i]) / dx - slope(data)[i])
        ds <- slope(data)[i+1] - slope(data)[i]
        secondDerivative <- (a - ds + ifelse(a>ds,1,-1)*sqrt((a-ds)^2+ds^2)) / dx
        breakingPoint <- (point.x(data)[i+1] + point.x(data)[i] + ds / secondDerivative) / 2

        leftBound[2*i-1] <- point.x(data)[i]
        rightBound[2*i-1] <- breakingPoint
        polynomial[[2*i-1]] <- point.y(data)[i] + (slope(data)[i] * polynomial(c(-point.x(data)[i],1))) + ((secondDerivative/2) * polynomial(c(point.x(data)[i]^2,-2*point.x(data)[i],1)))

        leftBound[2*i] <- breakingPoint
        rightBound[2*i] <- point.x(data)[i+1]
        polynomial[[2*i]] <- point.y(data)[i+1] + (slope(data)[i+1] * polynomial(c(-point.x(data)[i+1],1))) - ((secondDerivative/2) * polynomial(c(point.x(data)[i+1]^2,-2*point.x(data)[i+1],1)))
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}