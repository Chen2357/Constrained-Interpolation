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

setMethod("[", "pointData", function(x,i,...) pointData(x=point.x(x)[i],y=point.y(x)[i]))

setGeneric("graph", function(object,...) standardGeneric("graph"))
setMethod("graph", "pointData",
    function(object, col="red") {
        points(object@x,object@y,col=col)
    }
)

# ANCHOR: slopedPointData
slopedPointData <- setClass("slopedPointData",
    contains = "pointData",
    slots = c(
        slope = "vector"
    )
)

setValidity("slopedPointData", function(object) {
    if (length(object@x) != length(object@y)) {
        "Numbers of x and y data are mismatched"
    } else if (length(object@x) != length(object@slope)) {
        "Numbers of point and slope data are mismatched"
    } else {
        TRUE
    }
})

setGeneric("slope", function(object) standardGeneric("slope"))
setMethod("slope", "slopedPointData", function(object) object@slope)
setGeneric("slope<-", function(object,value) standardGeneric("slope<-"))
setMethod("slope<-", "slopedPointData", function(object,value) {
    object@slope <- value
    validObject(object)
    object
})