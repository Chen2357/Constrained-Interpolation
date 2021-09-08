#' An S4 Class to Represent a Polynomial
#' 
#' @details
#' Use \code{coef(object)} to access the coefficients.
#' Use \code{degree(object)} to access the degree of the polynomial.
#' Use \code{predict(object, newdata)} to evaluate the polynomial at \code{newdata} (polynomial is allowed). Alternatively, use \code{as.function(x)} to turn the polynomial into a function.
#' Use \code{differentiate(x)} to take the derivative of the polynomial.
#' 
#' @slot coef A vector of coefficients from the lowest degree term to the highest degree term. For example, \code{coef[1]} correspond to the constant term in the polynomial.
polynomial <- setClass("polynomial",
    slots = c(
        coef = "vector"
    )
)

setValidity("polynomial", function(object) {
    for (i in seq_len(length(object))) {
        if (!is.numeric(object@coef[i]) || is.na(object@coef[i])) {
            return(paste("Polynomial coefficients must be numeric, found", paste(object@coef, collapse = ", ")))
        }
    }
    TRUE
})

# ANCHOR Accessor
setMethod("coef", "polynomial", function(object) object@coef)
setMethod("coef<-", "polynomial", function(object,value) {
    object@coef <- value
    validObject(object)
    object
})

# ANCHOR Operator
setMethod("+", signature(e1 = "polynomial", e2 = "polynomial"), function(e1, e2) {
    n <- max(degree(e1), degree(e2))
    sum <- polynomial(degree = n)
    sum@coef <- c(coef(e1), rep(0, n-degree(e1))) + c(coef(e2), rep(0, n-degree(e2)))
    return(sum)
})

setMethod("+", signature(e1 = "numeric", e2 = "polynomial"), function(e1, e2) {
    coef(e2)[1] <- e1 + coef(e2)[1]
    return(e2)
})

setMethod("+", signature(e1 = "polynomial", e2 = "numeric"), function(e1, e2) {
    coef(e1)[1] <- coef(e1)[1] + e2
    return(e1)
})

setMethod("-", signature(e1 = "polynomial", e2 = "polynomial"), function(e1, e2) {
    n <- max(degree(e1), degree(e2))
    result <- polynomial(degree = n)
    result@coef <- c(coef(e1), rep(0, n-degree(e1))) - c(coef(e2), rep(0, n-degree(e2)))
    return(result)
})

setMethod("-", signature(e1 = "numeric", e2 = "polynomial"), function(e1, e2) {
    coef(e2) <- -coef(e2)
    coef(e2)[1] <- e1 + coef(e2)[1]
    return(e2)
})

setMethod("-", signature(e1 = "polynomial", e2 = "numeric"), function(e1, e2) {
    coef(e1)[1] <- coef(e1)[1] - e2
    return(e1)
})

setMethod("*", signature(e1 = "polynomial", e2 = "polynomial"), function(e1, e2) {
    product <- polynomial(degree = (degree(e1) + degree(e2)))

    for (i in seq_len(length(e1))) {
        for (j in seq_len(length(e2))) {
            product@coef[i+j-1] <- coef(product)[i+j-1] + coef(e1)[i] * coef(e2)[j]
        }
    }
    return(product)
})

setMethod("*", signature(e1 = "polynomial", e2 = "numeric"), function(e1, e2) return(polynomial(coef(e1) * e2)))


setMethod("*", signature(e1 = "numeric", e2 = "polynomial"), function(e1, e2) return(polynomial(e1 * coef(e2))))

setMethod("^", signature(e1 = "polynomial", e2 = "numeric"), function(e1, e2) {
    if (e2%%1 != 0) {
        warning("exponent must be an integer, caught ", e2)
        return(NULL)
    }
    result <- polynomial(c(1))

    while (e2 > 0) {
        if (e2 %% 2 == 1) {
            result <- result * e1
            e2 <- (e2-1)/2
        } else {
            e2 <- e2/2
        }
        e1 <- e1 * e1
    }

    return(result)
})

# ANCHOR Methods
setMethod("degree", "polynomial", function(x) length(x@coef)-1)
setMethod("degree<-", "polynomial", function(x,value) {
    length(x@coef) <- value + 1
    validObject(x)
    x
})

setMethod("length", "polynomial", function(x) length(x@coef))
setMethod("length<-", "polynomial", function(x,value) {
    length(x@coef) <- value
    validObject(x)
    x
})

setMethod("initialize", "polynomial",
    function(.Object, coef, degree) {
        if (missing(coef)) {
            .Object@coef <- rep(0,ifelse(missing(degree), 1, degree+1))
        } else if (missing(degree)) {
            .Object@coef <- coef
        } else {
            .Object@coef <- c(coef, rep(0, degree + 1 - length(coef)))
        }

        validObject(.Object)
        return(.Object)
    }
)

setMethod("differentiate", "polynomial", function(x) polynomial(coef(x)[-1] * 1:degree(x)))

setMethod("predict", signature(object="polynomial"),
    function(object, newdata) {
        if (class(newdata) == "numeric") {
            result <- rep(0, length(newdata))
            for (i in seq_len(length(object))) {
                result <- result + object@coef[i] * newdata ^ (i-1)
            }
            return(result)
        } else if (class(newdata) == "polynomial") {
            result <- polynomial(c(0))
            p <- 1
            for (i in seq_len(length(object))) {
                result <- result + object@coef[i] * p
                p <- p * newdata
            }
            return(result)
        } else if (class(newdata) == "dual") {
            result <- dual(0, degree = degree(newdata), length=length(newdata))
            p <- 1
            for (i in seq_len(length(object))) {
                result <- result + object@coef[i] * p
                p <- p * newdata
            }
            return(result)
        } else {
            stop("Unsupprted newdata class in predict where object is polynomial")
        }
    }
)

setMethod("as.function", "polynomial", function(x) function(xx) predict(x,xx))

setMethod("as.character", "polynomial",
    function(x, xlab="x", digits = getOption("digits")) {
        eq <- ""
        for (i in seq_len(length(x))) {
            if (x@coef[i] != 0) {
                if (i==1) {
                    eq <- paste(eq, signif(x@coef[i], digits), sep = "")
                }
                else {
                    eq <- paste(eq,ifelse(eq!=""," + ", ""),signif(x@coef[i], digits),"*",xlab,ifelse(i==2,"",paste("^",i-1,sep="")), sep = "")
                }
            }
        }
        return(ifelse(eq=="","0",eq))
    }
)

setMethod("show", "polynomial",
    function(object) {
        print(noquote(as.character(object)))
    }
)

# ANCHOR Extrema Functions

#' Finding the x-values of Extrema of a Polynomial
#' 
#' @param poly A `polynomial` type.
#' @param tol Tolerance.
#' @return A vector containing the x-values of the extrema.
polynomial.extrema.x <- function(poly, tol = sqrt(.Machine$double.eps)) {
    r <- polyroot(coef(differentiate(poly)))
    r <- Re(r[abs(Im(r)) < tol])
    i <- 1
    while(i <= length(r)) {
        I <- (abs(r - r[i]) < tol)
        if (sum(I) %% 2 == 1) {
            r <- r[!I | (seq_len(length(r)) == i)]
            i <- i + 1
        } else {
            r <- r[!I]
        }
    }
    r <- sort(r)
    return(r)
}

#' Finding the y-values of Extrema of a Polynomial
#' 
#' @param poly A `polynomial` type.
#' @param tol Tolerance.
#' @return A vector containing the y-values of the extrema.
polynomial.extrema.y <- function(poly, tol = sqrt(.Machine$double.eps)) {
    predict(poly, polynomial.extrema.x(poly, tol))
}

#' Finding the Extrema of a Polynomial
#' 
#' @param poly A `polynomial` type.
#' @param tol Tolerance.
#' @return A `pointData` type containing the extrema points.
polynomial.extrema <- function(poly, tol = sqrt(.Machine$double.eps)) {
    x <- polynomial.extrema.x(poly, tol)
    y <- predict(poly, x)
    return(pointData(x, y))
}