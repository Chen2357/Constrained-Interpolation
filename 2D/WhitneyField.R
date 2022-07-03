whitneyField <- setClass(
    "whitneyField",
    slots = c(
        x = "numeric",
        y = "numeric",
        coef = "matrix"
    )
)

setValidity("whitneyField", function(object) {
    if (length(object@x) != length(object@y)) {
        return("Numbers of x and y data are mismatched")
    } else if (ncol(object@coef) != 3) {
        return("Number of columns of coef must be 3")
    } else if (nrow(object@coef) != length(object@x)) {
        return("Number of x data and number of columns of coef mismatched")
    }
    return(TRUE)
})

setMethod("initialize", "whitneyField",
    function(.Object, x = numeric(0), y = numeric(0), coef = matrix(nrow=length(x), ncol=3)) {
        .Object@x <- x
        .Object@y <- y

        if (!is.matrix(coef)) {
            coef <- matrix(coef, nrow=length(x), ncol=3)
        }

        .Object@coef <- coef

        validObject(.Object)
        return(.Object)
    }
)