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
    function(.Object, x = numeric(0), y = numeric(0), w = numeric(0)) {
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
setMethod("[<-", "whitneySquare", function(x,i,...,value) {
    x@x[i] <- value@x
    x@y[i] <- value@y
    x@w[i] <- value@w
    return(x)
})

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

bisect.whitney <- function(square) {
    square@w <- c(square@w %*% t(c(1,1,1,1)/2))

    square@x <- c(square@x %*% t(c(1,1,1,1))) + square@w * c(0,1,0,1)
    square@y <- c(square@y %*% t(c(1,1,1,1))) + square@w * c(0,0,1,1)

    return(square)
}

whitneyDecomposition <- setClass("whitneyDecomposition",
    slots = c(
        squares = "whitneySquare",
        order = "numeric",
        root = "numeric",
        nodes = "matrix",
        children = "matrix"
    )
)
# setValidity("whitneyDecomposition", function(object) {
#     rownames(object@nodes) <- NULL
#     rownames(object@children) <- NULL

#     if (any(object@nodes[,1]==0)) {
#         print("Found error 1")
#         print(object)
#         return("Error 1")
#     }
#     print("Validity passed")
#     print(object)
#     return(TRUE)
# })

setMethod("initialize", "whitneyDecomposition",
    function(.Object) {
        .Object@squares <- whitneySquare()
        .Object@order <- 5
        .Object@root <- 0
        .Object@nodes <- matrix(nrow = 0, ncol = .Object@order-1)
        .Object@children <- matrix(nrow = 0, ncol = .Object@order)

        # validObject(.Object)
        return(.Object)
    }
)

insert.whitney <- function(decomposition, squares) {
    lastIndex <- length(decomposition@squares)
    decomposition@squares <- append(decomposition@squares, squares)

    order <- decomposition@order

    for (squareIndex in seq_len(length(squares))) {
        n <- decomposition@root
        key <- lastIndex + squareIndex

        if (n == 0) {
            decomposition@nodes <- rbind(decomposition@nodes, c(key, rep(0, order - 2)))
            decomposition@children <- rbind(decomposition@children, rep(0, order))
            decomposition@root <- nrow(decomposition@nodes)
            next
        }

        nodePath <- c()
        insertPath <- c()

        lastNotFull <- 0

        while (TRUE) {
            s <- decomposition@nodes[n,]
            i <- sum(decomposition@squares[s[s != 0]] < squares[squareIndex], na.rm=TRUE)
            notFull <- any(s == 0)

            nodePath <- c(nodePath, n)
            insertPath <- c(insertPath, i)
            
            if (notFull) lastNotFull <- length(nodePath)
            
            if (decomposition@children[n,1] == 0) {
                if (notFull) {
                    decomposition@nodes[n,] <- append(s[1:(order-2)], key, i)
                    lastNotFull <- NULL
                }
                break
            } else {
                n <- decomposition@children[n,i+1]
            }
        }

        if (!is.null(lastNotFull)) {
            medium <- ceiling(order/2)
            hold <- key
            holdChild <- 0

            for (level in length(nodePath):lastNotFull) {
                n <- nodePath[level]
                i <- insertPath[level]

                if (level != lastNotFull) {
                    extendedNode <- append(decomposition@nodes[n,], hold, i)
                    extendedChildren <- append(decomposition@children[n,], holdChild, i+1)

                    decomposition@nodes[n,] <- c(extendedNode[1:(medium-1)], rep(0, order-medium))
                    decomposition@children[n,] <- c(extendedChildren[1:medium], rep(0, order-medium))

                    newNode <- extendedNode[(medium+1):order]
                    newNode <- c(newNode, rep(0, medium - 1))
                    newChildren <- extendedChildren[(medium+1):(order+1)]
                    newChildren <- c(newChildren, rep(0, medium - 1))

                    decomposition@nodes <- rbind(decomposition@nodes, newNode)
                    decomposition@children <- rbind(decomposition@children, newChildren)

                    hold <- extendedNode[medium]
                    holdChild <- nrow(decomposition@nodes)
                } else if (level == 0) {
                    newNode <- c(hold, rep(0, order - 2))
                    newChildren <- c(decomposition@root, holdChild, rep(0, order - 2))

                    decomposition@nodes <- rbind(decomposition@nodes, newNode)
                    decomposition@children <- rbind(decomposition@children, newChildren)

                    decomposition@root <- nrow(decomposition@nodes)
                } else {
                    decomposition@nodes[n,] <- append(decomposition@nodes[n,-(order-1)], hold, i)
                    decomposition@children[n,] <- append(decomposition@children[n,-order], holdChild, i+1)
                }
            }
        }
    }

    rownames(decomposition@nodes) <- NULL
    rownames(decomposition@children) <- NULL

    return(decomposition)
}

partition.whitney <- function(x, y) {
    result <- insert.whitney(whitneyDecomposition(), whitneySquare(0,0,1))

    q <- 1
    queue <- c(1)

    while (q <= length(queue)) {
        # if (q > 65536) {
        #     print("Exceeds limit")
        #     break
        # }
        i <- queue[q]

        squarex <- result@squares@x[i]
        squarey <- result@squares@y[i]
        squarew <- result@squares@w[i]

        count <- sum((squarex - squarew <= x & x < squarex + 2 * squarew) & (squarey - squarew <= y & y < squarey + 2 * squarew))

        if (count > 1) {
            bisection <- bisect.whitney(result@squares[i])
            result@squares[i] <- bisection[1]

            result <- insert.whitney(result, bisection[2:4])
            queue <- c(queue, i, (length(result@squares)-2):length(result@squares))
        }
        q <- q + 1
    }

    return(result)
}

rect.whitney <- function(decomposition, x, y) {
    for (i in seq_len(length(decomposition@squares))) {
        square <- decomposition@squares[i]

        hasPoint <- any((square@x <= x & x < square@x + square@w) & (square@y <= y & y < square@y + square@w))
        col <- ifelse(hasPoint, "pink", "lightgreen")

        rect(square@x, square@y, square@x + square@w, square@y + square@w, col = col)
    }

    points(x, y, col = "red")
}

x <- c(0.2,0.2,0.4,0.3)
y <- c(0.3,0.1,0.5,0.9)

W <- partition.whitney(x, y)

plot(c(0,1), c(0,1), type = "n")
rect.whitney(W, x, y)