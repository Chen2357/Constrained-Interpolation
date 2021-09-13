#' An S4 Class to Represent a Piecewise Polynomial
#' 
#' @details
#' Use \code{length(object)} to get the number of piecewise ranges.
#' Use \code{predict(object, newdata)} to evaluate at \code{newdata}. Alternatively, use \code{as.function(x)} to turn the piecewise polynomial into a function.
#' Use \code{differentiate(x)} to take the derivative of each polynomial.
#' Use \code{plot(x)} to \code{lines{x}} to plot the piecewise polynomial.
#' 
#' @slot leftBound A vector that stores the left bounds of the piecewise ranges.
#' @slot rightBound A vector that stores the right bounds of the piecewise ranges.
#' @slot polynomial A list that stores the polynomials at each piecewise ranges.
piecewisePolynomial <- setClass("piecewisePolynomial",
    slots = c(
        leftBound = "vector",
        rightBound = "vector",
        polynomial = "list"
    )
)

setValidity("piecewisePolynomial", function(object) {
    if (length(object@leftBound) != length(object@rightBound)) {
        "Numbers of leftBound and rightBound are mismatched"
    } else if (length(object@leftBound) != length(object@polynomial)) {
        "Numbers of ranges and polynomials are mismatched"
    } else {
        for (i in seq_len(length(object@leftBound))) {
            if (object@leftBound[i] > object@rightBound[i]) {
                return("Left bound cannot be greater than right bound")
            }
            else if (object@leftBound[i] < object@leftBound && object@leftBound < object@rightBound[i]) {
                return("Ranges cannot overlap")
            }
        }
        TRUE
    }
})

## ANCHOR Operators

combinePiecewisePolynomial <- function(FUN, e1, e2, tol = sqrt(.Machine$double.eps)) {
    I1 <- sort.list(e1@leftBound)
    I2 <- sort.list(e2@leftBound)
    j1 <- 1
    j2 <- 1

    marker <- NULL
    resetMarker <- TRUE

    leftBound <- c()
    rightBound <- c()
    polynomial <- list()

    while((j1 <= length(I1)) & (j2 <= length(I2))) {
        i1 <- I1[j1]
        i2 <- I2[j2]

        if (resetMarker) {
            marker <- max(min(e1@leftBound[i1], e2@leftBound[i2]), marker)
            isMarker1 <- (e1@leftBound[i1] < e2@leftBound[i2])
            resetMarker <- FALSE
        }

        if (isMarker1) {
            leftBound2 <- e2@leftBound[i2]
            rightBound1 <- e1@rightBound[i1]
            rightBound2 <- e2@rightBound[i2]
            poly1 <- e1@polynomial[[i1]]
            poly2 <- e2@polynomial[[i2]]
        } else {
            leftBound2 <- e1@leftBound[i1]
            rightBound1 <- e2@rightBound[i2]
            rightBound2 <- e1@rightBound[i1]
            poly1 <- e2@polynomial[[i2]]
            poly2 <- e1@polynomial[[i1]]
        }
        
        if (rightBound1 < leftBound2 + tol) {
            leftBound <- c(leftBound, marker)
            rightBound <- c(rightBound, rightBound1)
            polynomial <- c(polynomial, poly1)

            if (isMarker1) j1 <- j1 + 1
            else j2 <- j2 + 1
            resetMarker <- TRUE
        } else {
            if (isMarker1) {
                result <- FUN(poly1, poly2)
            } else {
                result <- FUN(poly2, poly1)
            }

            if (abs(rightBound1 - rightBound2) < tol) {
                j1 <- j1 + 1
                j2 <- j2 + 1
                newMarker <- rightBound1
                resetMarker <- TRUE
            } else if (rightBound1 < rightBound2) {
                if (isMarker1) j1 <- j1 + 1
                else j2 <- j2 + 1
                newMarker <- rightBound1
                isMarker1 <- !isMarker1
            } else {
                if (!isMarker1) j1 <- j1 + 1
                else j2 <- j2 + 1
                newMarker <- rightBound2
            }

            if (abs(marker - leftBound2) > tol) {
                leftBound <- c(leftBound, marker, leftBound2)
                rightBound <- c(rightBound, leftBound2, newMarker)
                polynomial <- c(polynomial, poly1, result)
            } else {
                leftBound <- c(leftBound, leftBound2)
                rightBound <- c(rightBound, newMarker)
                polynomial <- c(polynomial, result)
            }

            marker <- newMarker
        }
    }

    if (!resetMarker) {
        leftBound <- c(leftBound, marker)
        if (isMarker1) {
            rightBound <- c(rightBound, e1@rightBound[I1[j1]])
            polynomial <- c(polynomial, e1@polynomial[I1[j1]])
            j1 <- j1 + 1
        } else {
            rightBound <- c(rightBound, e2@rightBound[I2[j2]])
            polynomial <- c(polynomial, e2@polynomial[I2[j2]])
            j2 <- j2 + 1
        }
    }

    if (j1 <= length(I1)) {
        I <- I1[j1:length(I1)]
        leftBound <- c(leftBound, e1@leftBound[I])
        rightBound <- c(rightBound, e1@rightBound[I])
        polynomial <- c(polynomial, e1@polynomial[I])
    } else if (j2 <= length(I2)) {
        I <- I2[j2:length(I2)]
        leftBound <- c(leftBound, e2@leftBound[I])
        rightBound <- c(rightBound, e2@rightBound[I])
        polynomial <- c(polynomial, e2@polynomial[I])
    }

    return(piecewisePolynomial(leftBound, rightBound, polynomial))
}

### ANCHOR Addition
setMethod("+", signature(e1 = "piecewisePolynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    combinePiecewisePolynomial(`+`, e1, e2)
})
setMethod("+", signature(e1 = "piecewisePolynomial", e2 = "polynomial"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x + e2)
    return(e1)
})
setMethod("+", signature(e1 = "piecewisePolynomial", e2 = "numeric"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x + e2)
    return(e1)
})
setMethod("+", signature(e1 = "polynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    e2@polynomial <- lapply(e2@polynomial, function(x) e1 + x)
    return(e2)
})
setMethod("+", signature(e1 = "numeric", e2 = "piecewisePolynomial"), function(e1, e2) {
    e2@polynomial <- lapply(e2@polynomial, function(x) e1 + x)
    return(e2)
})

### ANCHOR Multiplication
setMethod("*", signature(e1 = "piecewisePolynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    combinePiecewisePolynomial(`*`, e1, e2)
})
setMethod("*", signature(e1 = "piecewisePolynomial", e2 = "polynomial"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x * e2)
    return(e1)
})
setMethod("*", signature(e1 = "piecewisePolynomial", e2 = "numeric"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x * e2)
    return(e1)
})
setMethod("*", signature(e1 = "polynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    e2@polynomial <- lapply(e2@polynomial, function(x) e1 * x)
    return(e2)
})
setMethod("*", signature(e1 = "numeric", e2 = "piecewisePolynomial"), function(e1, e2) {
    e2@polynomial <- lapply(e2@polynomial, function(x) e1 * x)
    return(e2)
})

### ANCHOR Subtraction
setMethod("-", signature(e1 = "piecewisePolynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    combinePiecewisePolynomial(`-`, e1, e2)
})
setMethod("-", signature(e1 = "piecewisePolynomial", e2 = "polynomial"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x - e2)
    return(e1)
})
setMethod("-", signature(e1 = "piecewisePolynomial", e2 = "numeric"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x - e2)
    return(e1)
})
setMethod("-", signature(e1 = "polynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    e2@polynomial <- lapply(e2@polynomial, function(x) e1 - x)
    return(e2)
})
setMethod("-", signature(e1 = "numeric", e2 = "piecewisePolynomial"), function(e1, e2) {
    e2@polynomial <- lapply(e2@polynomial, function(x) e1 - x)
    return(e2)
})

### ANCHOR Division
setMethod("/", signature(e1 = "piecewisePolynomial", e2 = "numeric"), function(e1, e2) {
    e1@polynomial <- lapply(e1@polynomial, function(x) x / e2)
    return(e1)
})

### ANCHOR Append
setMethod("%+%", signature(e1 = "piecewisePolynomial", e2 = "piecewisePolynomial"), function(e1, e2) {
    result <- piecewisePolynomial(c(e1@leftBound, e2@leftBound), c(e1@rightBound, e2@rightBound), c(e1@polynomial, e2@polynomial))
    validObject(result)
    return(result)
})

## ANCHOR Generics

setMethod("length", "piecewisePolynomial", function(x) length(x@leftBound))

setMethod("predict", signature(object="piecewisePolynomial"),
    function(object,newdata) {
        if (class(newdata) == "numeric") {
            y <- rep(NA, length(newdata))
            for(i in seq_len(length(object))) {
                indices <- which(object@leftBound[i] <= newdata & newdata <= object@rightBound[i])
                y[indices] <- predict(object@polynomial[[i]],newdata[indices])
            }
            return(y)
        } else if (class(newdata) == "dual") {
            y <- dual(0, degree = degree(newdata), length=length(newdata))
            for(i in seq_len(length(object))) {
                I <- which(object@leftBound[i] <= newdata & newdata <= object@rightBound[i])
                y[I] <- predict(object@polynomial[[i]],newdata[I])
            }
            return(y)
        } else {
            stop("Unsupprted newdata class in predict where object is piecewisePolynomial")
        }
    }
)

setMethod("as.function", "piecewisePolynomial", function(x) function(xx) predict(x,xx))

setMethod("leftMostBound", "piecewisePolynomial", function(object) min(object@leftBound))

setMethod("rightMostBound", "piecewisePolynomial", function(object) max(object@rightBound))

setMethod("differentiate", "piecewisePolynomial", function(x) piecewisePolynomial(x@leftBound, x@rightBound, lapply(x@polynomial, differentiate)))

setMethod("plot", "piecewisePolynomial",
    function(x, interval=seq(leftMostBound(x),rightMostBound(x),0.05),type="l", ...) {
        plot(interval, predict(x, interval), type=type, ...)
    }
)

setMethod("lines", "piecewisePolynomial",
    function(x, interval=seq(leftMostBound(x),rightMostBound(x),0.05), ...) {
        lines(interval, predict(x, interval), ...)
    }
)

setMethod("initialize", "piecewisePolynomial",
    function(.Object, leftBound = numeric(0), rightBound = numeric(0), polynomial = list()) {
        .Object@leftBound <- leftBound
        .Object@rightBound <- rightBound
        .Object@polynomial <- polynomial

        validObject(.Object)
        return(.Object)
    }
)

defaultRangeFormatter <- function(min, max, xlab="x", digits = getOption("digits")) paste("(",xlab,">",signif(min,digits)," & ",xlab,"<",signif(max,digits),")", sep = "")

setMethod("as.character", "piecewisePolynomial",
    function(x, xlab="x", rangeFormatter = defaultRangeFormatter, digits = getOption("digits")) {
        eq <- ""
        for (i in seq_len(length(x))) {
            eq <- paste(eq,ifelse(eq=="",""," + "),rangeFormatter(x@leftBound[i],x@rightBound[i],xlab=xlab,digits=digits),"*(",as.character(x@polynomial[[i]],xlab=xlab,digits=digits),")", sep = "")
        }
        return(ifelse(eq=="","0",eq))
    }
)

setMethod("degree", "piecewisePolynomial", function(x) max(unlist(lapply(x@polynomial, degree), use.names=FALSE)))

setMethod("as.data.frame", "piecewisePolynomial", function(x, xlab="x", rangeFormatter = defaultRangeFormatter, digits = getOption("digits"), ...) {
    interval <- NULL
    equation <- NULL
    for (i in seq_len(length(x))) {
        interval <- append(interval, rangeFormatter(x@leftBound[i],x@rightBound[i],xlab=xlab, digits = digits))
        equation <- append(equation, as.character(x@polynomial[[i]],xlab=xlab,digits = digits))
    }
    return(data.frame(interval, equation, ...))
})

setMethod("show", "piecewisePolynomial",
    function(object) {
        print(noquote(as.character(object)))
    }
)

setMethod("as.piecewisePolynomial", "polynomial", function(object, leftBound, rightBound) {
    if (missing(leftBound)) leftBound <- -Inf
    if (missing(rightBound)) rightBound <- Inf
    return(piecewisePolynomial(leftBound, rightBound, list(object)))
})

setMethod("as.piecewisePolynomial", "numeric", function(object, leftBound, rightBound) {
    if (missing(leftBound)) leftBound <- -Inf
    if (missing(rightBound)) rightBound <- Inf
    piecewisePolynomial(leftBound, rightBound, list(polynomial(c(object))))
})

setMethod("as.piecewisePolynomial", "piecewisePolynomial", function(object, leftBound, rightBound) {
    if (missing(leftBound)) leftBound <- leftMostBound(object)
    if (missing(rightBound)) rightBound <- rightMostBound(object)
    leftBounds <- c()
    rightBounds <- c()
    polynomial <- list()
    
    for (i in seq_len(length(object))) {
        if (object@leftBound[i] < rightBound & object@rightBound[i] > leftBound) {
            leftBounds <- c(leftBounds, max(object@leftBound[i], leftBound))
            rightBounds <- c(rightBounds, min(object@rightBound[i], rightBound))
            polynomial <- c(polynomial, object@polynomial[[i]])
        }
    }

    return(piecewisePolynomial(leftBounds, rightBounds, polynomial))
})

# ANCHOR Extrema Functions

#' Finding the x-values of Extrema of a Piecewise Polynomial
#' 
#' @param poly A `piecewisePolynomial` type.
#' @param tol Tolerance.
#' @return A vector containing the x-values of the extrema.
piecewisePolynomial.extrema.x <- function(poly, tol = sqrt(.Machine$double.eps)) {
    point.x(piecewisePolynomial.extrema(poly, tol))
}

#' Finding the y-values of Extrema of a Piecewise Polynomial
#' 
#' @param poly A `piecewisePolynomial` type.
#' @param tol Tolerance.
#' @return A vector containing the y-values of the extrema.
polynomial.extrema.y <- function(poly, tol = sqrt(.Machine$double.eps)) {
    point.y(piecewisePolynomial.extrema(poly, tol))
}

#' Finding the Extrema of a Piecewise Polynomial
#' 
#' The ranges in this Piecewise Polynomial must be ordered.
#' The boundary points are also includede in the result.
#' 
#' @param poly A `piecewisePolynomial` type.
#' @param tol Tolerance.
#' @return A `pointData` type containing the extrema points.
piecewisePolynomial.extrema <- function(poly, tol = sqrt(.Machine$double.eps)) {
    x <- poly@leftBound[1]
    y <- predict(poly@polynomial[[1]], poly@leftBound[1])

    for (i in seq_len(length(poly))) {
        extrema <- polynomial.extrema(poly@polynomial[[i]])
        I <- ((point.x(extrema) > poly@leftBound[i] + tol) & (point.x(extrema) < poly@rightBound[i] - tol))
        x <- c(x, point.x(extrema)[I])
        y <- c(y, point.y(extrema)[I])

        if (i < length(poly)) {
            x1 <- poly@rightBound[i]
            x2 <- poly@leftBound[i+1]
            y1 <- predict(poly@polynomial[[i]], poly@rightBound[i])
            y2 <- predict(poly@polynomial[[i+1]], poly@leftBound[i+1])
            if ((x2 - x1 < tol) & (abs(y2 - y1) < tol)) {
                x <- c(x, x1)
                y <- c(y, y1)
            } else {
                x <- c(x, x1, x2)
                y <- c(y, y1, y2)
            }
        } else {
            x <- c(x, poly@rightBound[i])
            y <- c(y, predict(poly@polynomial[[i]], poly@rightBound[i]))
        }
    }

    return(pointData(x, y))
}