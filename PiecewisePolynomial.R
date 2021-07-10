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
        for(i in 1:length(object@leftBound)) {
            if(object@leftBound[i] > object@rightBound[i]) {
                return("Left bound cannot be greater than right bound")
            }
            else if (object@leftBound[i] < object@leftBound && object@leftBound < object@rightBound[i]) {
                return("Ranges cannot overlap")
            }
        }
        TRUE
    }
})

setMethod("length", "piecewisePolynomial", function(x) length(x@leftBound))

setMethod("func", "piecewisePolynomial",
    function(object) {
        function(x) {
            y = rep(NA, length(x))
            for(i in 1:length(object)) {
                indices <- which(object@leftBound[i] <= x & x <= object@rightBound[i])
                y[indices] = func(object@polynomial[[i]])(x[indices])
            }
            return(y)
        }
    }
)

setGeneric("leftMostBound", function(object) standardGeneric("leftMostBound"))
setMethod("leftMostBound", "piecewisePolynomial", function(object) min(object@leftBound))

setGeneric("rightMostBound", function(object) standardGeneric("rightMostBound"))
setMethod("rightMostBound", "piecewisePolynomial", function(object) max(object@rightBound))

setMethod("graph", "piecewisePolynomial",
    function(object,x=seq(leftMostBound(object),rightMostBound(object),0.05)) {
        plot(x, func(object)(x),type="l",xlab="x", ylab="y")
    }
)

setMethod("initialize", "piecewisePolynomial",
    function(.Object, leftBound, rightBound, polynomial) {
        .Object@leftBound <- leftBound
        .Object@rightBound <- rightBound
        .Object@polynomial <- polynomial

        validObject(.Object)
        return(.Object)
    }
)

defaultRangeFormatter <- function(min, max, x="x", digits = NULL) {
    if (!is.null(digits)) {
        min = format(min, digits=digits)
        max = format(max, digits=digits)
    }
    return(paste("(",x,">",min," & ",x,"<",max,")", sep = ""))
}

defaultPiecewisePolynomialFormat <- "raw"

setMethod("str", "piecewisePolynomial",
    function(object, x="x", rangeFormatter = defaultRangeFormatter, format = defaultPiecewisePolynomialFormat, digits = NULL) {
        if (format == "raw") {
            eq <- ""
            for (i in 1:length(object)) {
                eq <- paste(eq,ifelse(eq=="",""," + "),rangeFormatter(object@leftBound[i],object@rightBound[i],x=x,digits=digits),"*(",str(object@polynomial[[i]],x=x,digits=digits),")", sep = "")
            }
            return(ifelse(eq=="","0",eq))
        } else if (format == "table") {
            interval <- NULL
            equation <- NULL
            for (i in 1:length(object)) {
                interval <- append(interval, rangeFormatter(object@leftBound[i],object@rightBound[i],x=x,digits=digits))
                equation <- append(equation, str(object@polynomial[[i]],x=x,digits=digits))
            }
            result <- data.frame(interval, equation)
            return(result)
        } else {
            warning("Unknown defaultPiecewisePolynomialFormat: ", defaultPiecewisePolynomialFormat)
        }
    }
)

setMethod("show", "piecewisePolynomial",
    function(object) {
        print(str(object))
    }
)

setMethod("print", "piecewisePolynomial",
    function(x, ...) {
        print(str(x, ...))
    }
)