# A whitneySquare instance represents a collection of squares in a plane
# 
# `x` is a vector of the x-coordinates of the lower left corners of the squares.
# `y` is a vector of the y-coordinates of the lower left corners of the squares.
# `w` is a vector of the width of the squares.
# `field` is a matrix whose row n represents the linear polynomial associated with the square at index n.
whitneySquare <- setClass("whitneySquare",
    slots = c(
        x = "numeric",
        y = "numeric",
        w = "numeric",
        field = "matrix"
    )
)

setValidity("whitneySquare", function(object) {
    if (length(object@x) != length(object@y)) {
        return("Numbers of x and y data are mismatched")
    } else if (length(object@x) != length(object@w)) {
        return("Numbers of x and w data are mismatched")
    } else if (ncol(object@field) != 3) {
        return("Number of columns of field must be 3")
    } else if (nrow(object@field) != length(object@x)) {
        return("Number of x data and number of columns of field mismatched")
    }
    return(TRUE)
})

setMethod("initialize", "whitneySquare",
    function(.Object, x = numeric(0), y = numeric(0), w = numeric(0), field = matrix(nrow=length(x), ncol=3)) {
        .Object@x <- x
        .Object@y <- y
        .Object@w <- w
        .Object@field <- field

        validObject(.Object)
        return(.Object)
    }
)

setMethod("length", "whitneySquare", function(x) length(x@x))

setMethod("append", signature(x = "whitneySquare", values = "whitneySquare"), function(x, values, after = length(x)) {
    return(whitneySquare(append(x@x, values@x, after), append(x@y, values@y, after), append(x@w, values@w, after), rbind(x@field, values@field)))
})

setMethod("[", "whitneySquare", function(x,i,...) whitneySquare(x@x[i], x@y[i], x@w[i], x@field[i,]))
setMethod("[<-", "whitneySquare", function(x,i,...,value) {
    x@x[i] <- value@x
    x@y[i] <- value@y
    x@w[i] <- value@w
    x@field[i,] <- value@field[i,]
    return(x)
})

# `bisect.whitney` partitions each of the squares in the whitneySquare instance into four equal sized whitneySquare.
# The partition happens in-place, meaning the four squares resulted from one parent square will have adjacent index.
# It returns another whitneySquare instance with four times the number of squares.
bisect.whitney <- function(square) {
    square@w <- c(c(0.5,0.5,0.5,0.5) %*% t(square@w))

    square@x <- c(c(1,1,1,1) %*% t(square@x)) + square@w * c(0,1,0,1)
    square@y <- c(c(1,1,1,1) %*% t(square@y)) + square@w * c(0,0,1,1)

    square@field <- square@field[rep(seq_len(square), each = 2), ]

    return(square)
}