#' An S4 Class to Represent Patching Function
#' 
#' @details
#' A well behaved patching function `theta` should evaluate to 1 at 0, and 0 at 1.
#' 
#' @slot theta A function in the format of `function(x)` that returns a numeric value.
#' @slot description A function in the format of `function(a, b)` where `a` and `b` are left and right bounds of an interval and returns a string.
patchingFunction <- setClass(
    "patchingFunction",
    slots = c(
        theta = "function"
    )
)

setMethod("initialize", "patchingFunction",
    function(.Object, theta) {
        .Object@theta <- theta
        return(.Object)
    }
)

#' Polynomial Patching Function
#' 
#' @details
#' It inherits from `patchingFunction` class.
#' 
#' @slot thetaPolynomial A `polynomial` type.
patchingPolynomial <- setClass(
    "patchingPolynomial",
    slots = c(
        polynomial = "polynomial"
    ),
    contains = "patchingFunction"
)

setMethod("initialize", "patchingPolynomial",
    function(.Object, polynomial) {
        .Object@polynomial <- polynomial
        .Object@theta <- as.function(polynomial)
        return(.Object)
    }
)

patch.linear <- patchingPolynomial(polynomial(c(1, -1)))
patch.cubic <- patchingPolynomial(polynomial(c(1, 0, -3, 2)))
patch.fifthDegree <- patchingPolynomial(polynomial(c(1, 0, 0,-10, 15, -6)))
patch.bump <- patchingFunction(function(x) exp(1-1/(1-x^2)))

patching <- setClass(
    "patching",
    slots = c(
        func = "list",
        breaks = "vector",
        patch = "patchingFunction"
    )
)

setMethod("initialize", "patching",
    function(.Object, func = list(), breaks = numeric(0), patch) {
        .Object@func <- func
        .Object@breaks <- breaks
        .Object@patch <- patch

        validObject(.Object)
        return(.Object)
    }
)

setMethod("predict", signature(object="patching"),
    function(object, newdata) {
        n <- length(object@breaks)

        if (class(newdata) == "numeric") {
            result <- numeric(length(newdata))
        } else if (class(newdata) == "dual") {
            result <- dual(0, degree = degree(newdata), length = length(newdata))
        }

        if (n == 1) return(predict(object@func[[1]], newdata))

        I <- which(newdata <= object@breaks[1])
        result[I] <- predict(object@func[[1]], newdata[I])

        I <- which(newdata > object@breaks[n])
        result[I] <- predict(object@func[[n]], newdata[I])

        for (i in seq_len(n-1)) {
            I <- which(object@breaks[i] < newdata & newdata <= object@breaks[i+1])

            p <- (newdata[I] - object@breaks[i]) / (object@breaks[i+1] - object@breaks[i])
            a <- predict(object@func[[i]], newdata[I])
            b <- predict(object@func[[i+1]], newdata[I])

            result[I] <- object@patch@theta(p) * (a - b) + b
        }
        return(result)
    }
)

setMethod("as.piecewisePolynomial", "patching", function(object, leftBound, rightBound, tol = sqrt(.Machine$double.eps)) {
    if (class(object@patch) != "patchingPolynomial") return(NULL)

    if (missing(leftBound)) leftBound <- -Inf
    if (missing(rightBound)) rightBound <- Inf
    n <- length(object@breaks)
    result <- piecewisePolynomial()

    if (n==1) return(as.piecewisePolynomial(object@func[[1]], leftBound, rightBound))

    if (leftBound < object@breaks[1]) {
        r <- min(object@breaks[1], rightBound)
        poly <- as.piecewisePolynomial(object@func[[1]], leftBound, r)
        if (is.null(poly)) return(NULL)
        result <- result %+% poly
    }

    for (i in seq_len(n-1)) {
        if (leftBound < object@breaks[i+1] - tol & rightBound > object@breaks[i] + tol) {
            l <- max(leftBound, object@breaks[i])
            r <- min(rightBound, object@breaks[i+1])
        } else if (rightBound < object@breaks[i]) {
            return(result)
        } else {
            next
        }

        a <- as.piecewisePolynomial(object@func[[i]], l, r)
        b <- as.piecewisePolynomial(object@func[[i+1]], l, r)
        p <- percentagePolynomial(object@breaks[i], object@breaks[i+1])
        if (is.null(a) | is.null(b)) return(NULL)

        result <- result %+% (object@patch@theta(p) * (a - b) + b)
    }

    if (object@breaks[n] < rightBound) {
        l <- max(object@breaks[n], leftBound)
        poly <- as.piecewisePolynomial(object@func[[n]], l, rightBound)
        if (is.null(poly)) return(NULL)
        result <- result %+% poly
    }

    return(result)
})

setMethod("+", signature(e1 = "patching", e2 = "numeric"), function(e1, e2) {
    e1@func <- lapply(e1@func, function(x) x + e2)
    return(e1)
})
setMethod("+", signature(e1 = "numeric", e2 = "patching"), function(e1, e2) {
    e2@func <- lapply(e2@func, function(x) e1 + x)
    return(e2)
})

setMethod("*", signature(e1 = "patching", e2 = "numeric"), function(e1, e2) {
    e1@func <- lapply(e1@func, function(x) x * e2)
    return(e1)
})
setMethod("*", signature(e1 = "numeric", e2 = "patching"), function(e1, e2) {
    e2@func <- lapply(e2@func, function(x) e1 * x)
    return(e2)
})

setMethod("/", signature(e1 = "patching", e2 = "numeric"), function(e1, e2) {
    e1@func <- lapply(e1@func, function(x) x / e2)
    return(e1)
})