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
        theta = "function",
        description = "function"
    )
)

setMethod("initialize", "patchingFunction",
    function(.Object, theta, description) {
        .Object@theta <- theta
        .Object@description <- description

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
        .Object@description <- function(a, b) as.character(predict(polynomial, percentagePolynomial(a, b)))

        return(.Object)
    }
)

patchingDifferentiable <- setClass(
    "patchingDifferentiable",
    slots = c(
        derivative = "patchingFunction"
    ),
    contains = "patchingFunction"
)

cubicPatch <- patchingPolynomial(polynomial(c(1, 0, -3, 2)))
fifthDegreePatch <- patchingPolynomial(polynomial(c(1, 0, 0,-10, 15, -6)))
bumpPatch <- patchingFunction(
    theta = function(x) exp(1-1/(1-x^2)),
    description = function(a, b) paste("e^(1-1/(", 1-(polynomial(c(-a/(b-a), 1/(b-a)))^2), "))", sep="")
)

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

        I <- which(newdata <= object@breaks[1])
        result[I] <- predict(object@func[[1]], newdata[I])

        I <- which(newdata > object@breaks[n])
        result[I] <- predict(object@func[[n]], newdata[I])

        if (n == 1) return(result)

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
