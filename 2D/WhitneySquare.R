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

        if (!is.matrix(field)) {
            field <- matrix(field, nrow=length(x), ncol=3)
        }

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

setMethod("is.na", signature(x = "whitneySquare"), function(x) {
    return(is.na(x@x) | is.na(x@y) | is.na(x@w))
})

# `bisect.whitney` partitions each of the squares in the whitneySquare instance into four equal sized whitneySquare.
# The partition happens in-place, meaning the four squares resulted from one parent square will have adjacent index.
# It returns another whitneySquare instance with four times the number of squares.
bisect.whitney <- function(square) {
    square@field <- square@field[rep(seq_len(length(square)), each = 4), ]

    square@w <- c(c(0.5,0.5,0.5,0.5) %*% t(square@w))

    square@x <- c(c(1,1,1,1) %*% t(square@x)) + square@w * c(0,1,0,1)
    square@y <- c(c(1,1,1,1) %*% t(square@y)) + square@w * c(0,0,1,1)

    return(square)
}

search.whitney <- function(squares, x, y, na.rm = FALSE) {
    result <- whitneySquare()
    for (i in seq_along(x)) {
        contain <- (squares@x <= x[i] & x[i] < squares@x + squares@w) & (squares@y <= y[i] & y[i] < squares@y + squares@w)
        n <- which(contain)[1]
        if (is.na(n)) {
            if (!na.rm) result <- append(result, whitneySquare(NA_integer_, NA_integer_, NA_integer_))
        }
        else {
            result <- append(result, squares[n])
        }
    }
    return(result)
}

partition.whitney <- function(field) {
    squares <- whitneySquare()
    queue <- whitneySquare(0,0,1)

    x <- field@x
    y <- field@y

    while (length(queue) > 0) {
        # TODO This part need to be modified for readability
        contain_matrix <- outer(queue@x - queue@w, field@x, `<=`) & outer(queue@x + 2 * queue@w, field@x, `>`) & outer(queue@y - queue@w, field@y, `<=`) & outer(queue@y + 2 * queue@w, field@y, `>`)

        contain_indices <- apply(contain_matrix, 1, function(x) which(x)[1])

        contain_count <- rowSums(contain_matrix)
        contain_zero <- contain_count == 0
        contain_one <- contain_count == 1

        ok <- contain_one | contain_zero

        queue@field[contain_one,] <- field@coef[contain_indices[contain_one],]

        queue@field[!ok,] <- matrix(field@coef[contain_indices[!ok],], ncol = 3)

        squares <- append(squares, queue[ok])
        queue <- bisect.whitney(queue[!ok])
    }

    return(squares)
}

rect.whitney <- function(squares, x, y) {
    for (i in seq_along(squares)) {
        square <- squares[i]

        hasPoint <- any((square@x <= x & x < square@x + square@w) & (square@y <= y & y < square@y + square@w))
        col <- ifelse(hasPoint, "pink", "lightgreen")

        rect(square@x, square@y, square@x + square@w, square@y + square@w, col = col)

        text(square@x + square@w/2, square@y + square@w/2, labels = square@field[1,1])
    }

    points(x, y, col = "red")
}