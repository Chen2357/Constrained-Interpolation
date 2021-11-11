whitneySquare <- setClass("whitneySquare",
    slots = c(
        x = "numeric",
        y = "numeric",
        w = "numeric"
    )
)

setValidity("whitneySquare", function(object) {
    if (length(object@x) != length(object@y)) {
        return("Numbers of x and y data are mismatched")
    } else if (length(object@x) != length(object@w)) {
        return("Numbers of x and w data are mismatched")
    }
    return(TRUE)
})

setMethod("initialize", "whitneySquare",
    function(.Object, x, y, w) {
        .Object@x <- x
        .Object@y <- y
        .Object@w <- w

        validObject(.Object)
        return(.Object)
    }
)

setMethod("length", "whitneySquare", function(x) length(x@x))

setMethod("append", signature(x = "whitneySquare", values = "whitneySquare"), function(x, values, after = length(x)) {
    return(whitneySquare(append(x@x, values@x, after), append(x@y, values@y, after), append(x@w, values@w, after)))
})

setMethod("[", "whitneySquare", function(x,i,...) whitneySquare(x@x[i], x@y[i], x@w[i]))

setMethod(">", signature(e1 = "whitneySquare", e2 = "whitneySquare"), function(e1, e2) {
    result <- rep(FALSE, max(length(e1), length(e2)))

    result[e1@y > e2@y] <- TRUE
    result[e1@y < e2@y] <- FALSE

    I <- e1@y == e2@y
    result[I] <- (e1@x > e2@x)[I]

    return(result)
})

setMethod("<", signature(e1 = "whitneySquare", e2 = "whitneySquare"), function(e1, e2) {
    result <- rep(FALSE, max(length(e1), length(e2)))

    result[e1@y < e2@y] <- TRUE
    result[e1@y > e2@y] <- FALSE

    I <- e1@y == e2@y
    result[I] <- (e1@x < e2@x)[I]

    return(result)
})

whitneyDecomposition <- setClass("whitneyDecomposition",
    slots = c(
        squares = "whitneySquare",
        order = "numeric",
        root = "numeric",
        nodes = "matrix",
        children = "matrix"
    )
)


