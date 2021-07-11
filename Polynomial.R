polynomial <- setClass("polynomial",
    slots = c(
        coef = "vector"
    )
)

setValidity("polynomial", function(object) {
    for (i in 1:length(object@coef)) {
        if(!is.numeric(object@coef[i]) || is.na(object@coef[i])) {
            print(object@coef)
            return("Polynomial coefficients must be numeric")
        }
    }
    TRUE
})

# ANCHOR Accessor
setGeneric("coef", function(object) standardGeneric("coef"))
setMethod("coef", "polynomial", function(object) object@coef)
setGeneric("coef<-", function(object,value) standardGeneric("coef<-"))
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

    for (i in 1:length(coef(e1))) {
        for (j in 1:length(coef(e2))) {
            product@coef[i+j-1] <- coef(product)[i+j-1] + coef(e1)[i] * coef(e2)[j]
        }
    }
    return(product)
})

setMethod("*", signature(e1 = "polynomial", e2 = "numeric"), function(e1, e2) return(polynomial(coef(e1) * e2)))


setMethod("*", signature(e1 = "numeric", e2 = "polynomial"), function(e1, e2) return(polynomial(e1 * coef(e2))))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

setMethod("^", signature(e1 = "polynomial", e2 = "numeric"), function(e1, e2) {
    if(e2%%1 != 0) {
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

setGeneric("chain", function(e1,e2) standardGeneric("chain"))
setMethod("chain", signature(e1 = "polynomial", e2 = "polynomial"), function(e1, e2) {
    result <- polynomial(c(0))
    for (i in 1:length(coef(e1))) {
        result <- result + (polynomial(c(coef(e1)[i])) * (e2 ^ (i-1)))
    }
    return(result)
})

# ANCHOR Methods
setGeneric("degree", function(object) standardGeneric("degree"))
setMethod("degree", "polynomial", function(object) length(object@coef)-1)
setGeneric("degree<-", function(object,value) standardGeneric("degree<-"))
setMethod("degree<-", "polynomial", function(object,value) {
    length(object@coef) <- value + 1
    validObject(object)
    object
})

setMethod("initialize", "polynomial",
    function(.Object, coef, degree) {
        if(missing(coef)) {
            .Object@coef <- rep(0,ifelse(missing(degree), 1, degree+1))
        } else if(missing(degree)) {
            .Object@coef <- coef
        } else {
            .Object@coef <- c(coef, rep(0, degree + 1 - length(coef)))
        }

        validObject(.Object)
        return(.Object)
    }
)

setGeneric("func", function(object) standardGeneric("func"))
setMethod("func", "polynomial",
    function(object) {
        function(x) {
            result <- rep(0, length(x))
            for (i in 1:length(object@coef)) {
                result <- result + object@coef[i] * x ^ (i-1)
            }
            return(result)
        }
    }
)

setMethod("str", "polynomial",
    function(object, x="x", digits = NULL) {
        eq <- ""
        for (i in 1:length(object@coef)) {
            if (object@coef[i] != 0) {
                if (is.null(digits)) {
                    coef <- object@coef[i]
                } else {
                    coef <- format(object@coef[i], digits)
                }
                if (i==1) {
                    eq <- paste(eq, coef, sep = "")
                }
                else {
                    eq <- paste(eq,ifelse(eq!=""," + ", ""),coef,"*",x,ifelse(i==2,"",paste("^",i-1,sep="")), sep = "")
                }
            }
        }
        return(ifelse(eq=="","0",eq))
    }
)

setMethod("show", "polynomial",
    function(object) {
        print(str(object))
    }
)

# ANCHOR Make numeric conforms to polynomial
# setMethod("coef", "numeric", function(object) object)
# setMethod("degree", "numeric", function(object) 0)