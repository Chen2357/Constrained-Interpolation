dual <- setClass(
    "dual",
    slots = c(
        values = "vector",
        degree = "numeric"
    )
)

setValidity("dual", function(object) {
    if (length(object@degree) != 1) {
        return("Dual degree must be a single number")
    }
    if (object@degree < 0 || !is.wholenumber(object@degree)) {
        return("Dual degree must be a nonnegative integer")
    }
    if (length(object@values) %% (object@degree+1) != 0) {
        return("Dual values count must be a multiple of (degree+1)")
    }
    return(TRUE)
})

setMethod("initialize", "dual",
    function(.Object, values = numeric(), degree = length(values)-1, bydegree = FALSE) {
        remain <- -length(values) %% (degree+1)
        values <- c(values, rep(0, remain))
        if (bydegree) {
            values <- as.vector(t(matrix(values, ncol=degree+1)))
        }
        .Object@values <- values
        .Object@degree <- degree

        validObject(.Object)
        return(.Object)
    }
)

setMethod("length", "dual", function(x) length(x@values) / (x@degree+1))
setMethod("degree", "dual", function(x) x@degree)

setMethod("[", "dual", function(x,i,...) dual(x@values[i %x% rep(x@degree+1,x@degree+1) - x@degree:0], x@degree))
setMethod("[<-", "dual", function(x,i,...,value) {
    if (degree(value) != degree(x)) stop("Dual degrees do not match")
    x@values[i %x% rep(x@degree+1,x@degree+1) - x@degree:0] <- value@values
    validObject(x)
    x
})

setMethod("[[", "dual", function(x,i=seq_len(length(x)),j=0,...) {
    if (length(j) != 1) stop("Degree must be a single number")
    x@values[(i-1) * (x@degree+1) + j + 1]
})
setMethod("[[<-", "dual", function(x,i=seq_len(length(x)),j=0,...,value) {
    if (length(j) != 1) stop("Degree must be a single number")
    x@values[(i-1) * (x@degree+1) + j + 1] <- value
    validObject(x)
    x
})

# ANCHOR Addition
setMethod("+", signature(e1 = "dual", e2 = "dual"), function(e1, e2) {
    if (degree(e1) != degree(e2)) stop("Dual degrees do not match")
    e1@values <- e1@values + e2@values
    return(e1)
})
setMethod("+", signature(e1 = "dual", e2 = "numeric"), function(e1, e2) {
    if (length(e1) < length(e2)) stop("Length of the dual must be longer than length of vector")
    e1[[,0]] <- e1[[,0]] + e2
    return(e1)
})
setMethod("+", signature(e1 = "numeric", e2 = "dual"), function(e1, e2) {
    if (length(e1) > length(e2)) stop("Length of the dual must be longer than length of vector")
    e2[[,0]] <- e1 + e2[[,0]]
    return(e2)
})

## ANCHOR Multiplication
setMethod("*", signature(e1 = "dual", e2 = "dual"), function(e1, e2) {
    if (degree(e1) != degree(e2)) stop("Dual degrees do not match")
    n <- max(length(e1), length(e2))
    if (n %% length(e1) != n %% length(e2)) warning("longer object length is not a multiple of shorter object length")
    stride <- degree(e1) + 1
    values <- rep(0, n * stride)
    for (i in 0:degree(e1)) {
        I <- seq(i+1, (n-1) * stride + i + 1, stride)
        for (j in 0:i) {
            values[I] <- values[I] + e1[[,j]] * e2[[,i-j]]
        }
    }
    return(dual(values, degree(e1)))
})
setMethod("*", signature(e1 = "dual", e2 = "numeric"), function(e1, e2) {
    if (length(e1) < length(e2)) stop("Length of the dual must be longer than length of vector")
    e1@values <- e1@values * c(e2 %x% rep(1,degree(e1)+1))
    return(e1)
})
setMethod("*", signature(e1 = "numeric", e2 = "dual"), function(e1, e2) {
    if (length(e1) > length(e2)) stop("Length of the dual must be longer than length of vector")
    return(e2)
})

## ANCHOR Division
setMethod("/", signature(e1 = "dual", e2 = "numeric"), function(e1, e2) {
    if (length(e1) < length(e2)) stop("Length of the dual must be longer than length of vector")
    e1@values <- e1@values / c(e2 %x% rep(1,degree(e1)+1))
    return(e1)
})

setMethod("/", signature(e1 = "numeric", e2 = "dual"), function(e1, e2) {
    if (length(e1) > length(e2)) stop("Length of the dual must be longer than length of vector")
    n <- length(e2)
    if (n %% length(e1) != 0) warning("longer object length is not a multiple of shorter object length")
    stride <- degree(e2) + 1
    values <- rep(0, n * stride)
    I0 <- seq(1, (n-1) * stride + 1, stride)
    values[I0] <- e1 / e2@values[I0]

    for (i in 1:degree(e2)) {
        I <- seq(i+1, (n-1) * stride + i + 1, stride)
        sum <- rep(0, n)
        for (j in 0:(i-1)) {
            J <- seq(j+1, (n-1) * stride + j + 1, stride)
            sum <- sum + values[J] * e2[[,i-j]]
        }
        values[I] <- -sum / e2@values[I0]
    }
    return(dual(values, degree(e2)))
})

setMethod("exp", signature(x = "dual"), function(x) {
    n <- length(x)
    stride <- degree(x) + 1
    values <- rep(0, n * stride)
    I <- seq(1, (n-1) * stride + 1, stride)
    values[I] <- exp(x@values[I])

    if (degree(x) < 1) return(dual(values, degree(x)))

    binom <- rep(0, stride)
    binom[1] <- 1

    for (i in 1:degree(x)) {
        I <- seq(i+1, (n-1) * stride + i + 1, stride)
        f <- 1
        for (j in 0:(i-1)) {
            J <- seq(j+1, (n-1) * stride + j + 1, stride)
            values[I] <- values[I] + binom[j+1] * x[[,i-j]] * values[J] / f
            f <- f * (i-j)
        }

        binom <- binom + c(0, binom[1:(stride-1)])
    }

    return(dual(values, degree(x)))
})

setMethod("as.character", "dual",
    function(x, lab = "e", digits = getOption("digits")) {
        result <- rep("", length(x))
        for (i in 0:degree(x)) {
            I <- which(x[[,i]] != 0)
            nonempty <- I & which(result != "")
            if (i==0) {
                result[I] <- paste(result[I], signif(x[[I,i]], digits), sep = "")
            } else {
                result[nonempty] <- paste(result[nonempty], " + ", sep = "")
                result[I] <- paste(result[I], signif(x[[I,i]], digits),"*",lab,ifelse(i==1,"",paste("^",i,sep="")), sep = "")
            }
        }
        return(result)
    }
)

setMethod("show", "dual",
    function(object) {
        print(noquote(paste(as.character(object), collapse = ",  ")))
    }
)