# ANCHOR Data Point
setGeneric("point.x", function(object) standardGeneric("point.x"))
setGeneric("point.x<-", function(object,value) standardGeneric("point.x<-"))

setGeneric("point.y", function(object) standardGeneric("point.y"))
setGeneric("point.y<-", function(object,value) standardGeneric("point.y<-"))

# ANCHOR Matrix
setGeneric("symmetrize", function(x) standardGeneric("symmetrize"))

# ANCHOR Polynomial
setGeneric("coef", function(object) standardGeneric("coef"))
setGeneric("coef<-", function(object,value) standardGeneric("coef<-"))

setGeneric("degree", function(x) standardGeneric("degree"))
setGeneric("degree<-", function(x,value) standardGeneric("degree<-"))

# ANCHOR Piecewise functions
setGeneric("leftMostBound", function(object) standardGeneric("leftMostBound"))
setGeneric("rightMostBound", function(object) standardGeneric("rightMostBound"))

setGeneric("as.piecewisePolynomial", function(object, leftBound, rightBound) standardGeneric("as.piecewisePolynomial"))

# ANCHOR Differentiation
setGeneric("differentiate", function(x) standardGeneric("differentiate"))

## ANCHOR Interpolation

#' Test if a number is an integer
#' 
#' @param x A number
#' @param tol Tolerance
#' @return Whether or not `x` is an integer within tolerance
is.wholenumber <- function(x, tol = sqrt(.Machine$double.eps))  abs(x - round(x)) < tol