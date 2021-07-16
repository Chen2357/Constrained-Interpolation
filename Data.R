#' An S4 Class to Represent a Points
#' 
#' @details
#' Use \code{length(object)} to get the number of points.
#' Use \code{point.x(object)} or \code{point.y(object)} to access the x-coordinates or y-coordinates, respectively.
#' Use \code{differentiate(x)} to take the derivative of each polynomial.
#' Use \code{points(x)} to plot the points.
#' Use \code{[]} to access a subset of points.
#' 
#' @slot x A vector that stores x-coordinates.
#' @slot y A vector that stores y-coordinates.
pointData  <- setClass("pointData",
    slots = c(
        x = "vector",
        y = "vector"
    )
)

setValidity("pointData", function(object) {
    if (length(object@x) != length(object@y)) {
        "Numbers of x and y data are mismatched"
    } else if (is.unsorted(x)) {
        "x must be sorted"
    } else {
        TRUE
    }
})

setMethod("length", "pointData", function(x) length(x@x))

setGeneric("point.x", function(object) standardGeneric("point.x"))
setMethod("point.x", "pointData", function(object) object@x)
setGeneric("point.x<-", function(object,value) standardGeneric("point.x<-"))
setMethod("point.x<-", "pointData", function(object,value) {
    object@x <- value
    validObject(object)
    object
})

setGeneric("point.y", function(object) standardGeneric("point.y"))
setMethod("point.y", "pointData", function(object) object@y)
setGeneric("point.y<-", function(object,value) standardGeneric("point.y<-"))
setMethod("point.y<-", "pointData", function(object,value) {
    object@y <- value
    validObject(object)
    object
})

setMethod("points", "pointData",
    function(x, ...) {
        points(x@x,x@y,...)
    }
)

setMethod("[", "pointData", function(x,i,...) pointData(x=point.x(x)[i],y=point.y(x)[i]))
