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
`%+%` <- function(poly1, poly2) {
    n <- max(degree(poly1), degree(poly2))
    sum <- polynomial(degree = n)
    sum@coef <- c(coef(poly1), rep(0, n-degree(poly1))) + c(coef(poly2), rep(0, n-degree(poly2)))
    return(sum)
}

`%-%` <- function(poly1, poly2) {
    n <- max(degree(poly1), degree(poly2))
    result <- polynomial(degree = n)
    result@coef <- c(coef(poly1), rep(0, n-degree(poly1))) - c(coef(poly2), rep(0, n-degree(poly2)))
    return(result)
}

`%*%` <- function(poly1, poly2) {
    product <- polynomial(degree = (degree(poly1) + degree(poly2)))

    for (i in 1:length(coef(poly1))) {
        for (j in 1:length(coef(poly2))) {
            product@coef[i+j-1] <- coef(product)[i+j-1] + coef(poly1)[i] * coef(poly2)[j]
        }
    }
    return(product)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

`%^%` <- function(poly, exp) {
    if(exp%%1 != 0) {
        warning("exponent must be an integer, caught ", exp)
        return(NULL)
    }
    result <- polynomial(c(1))

    while (exp > 0) {
        if (exp %% 2 == 1) {
            result <- result %*% poly
            exp <- (exp-1)/2
        } else {
            exp <- exp/2
        }
        poly <- poly %*% poly
    }

    return(result)
}

`%chain%` <- function(poly1, poly2) {
    result <- polynomial(c(0))
    for (i in 1:length(coef(poly1))) {
        result <- result %+% (polynomial(c(coef(poly1)[i])) %*% (poly2 %^% (i-1)))
    }
    return(result)
}

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
                    coef <- object@coef[1]
                } else {
                    coef <- format(object@coef[1], digits)
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
setMethod("coef", "numeric", function(object) object)
setMethod("degree", "numeric", function(object) 0)