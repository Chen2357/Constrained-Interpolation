library(microbenchmark)

whitneySquare <- setClass("whitneySquare",
    slots = c(
        x = "numeric",
        y = "numeric",
        w = "numeric"
    )
)

setMethod("initialize", "whitneySquare",
    function(.Object, x = numeric(0), y = numeric(0), w = numeric(0)) {
        .Object@x <- x
        .Object@y <- y
        .Object@w <- w
        return(.Object)
    }
)

setMethod("[", "whitneySquare", function(x,i,...) whitneySquare(x@x[i], x@y[i], x@w[i]))
setMethod("[<-", "whitneySquare", function(x,i,...,value) {
    x@x[i] <- value@x
    x@y[i] <- value@y
    x@w[i] <- value@w
    return(x)
})

setMethod("length", "whitneySquare", function(x) length(x@x))

setMethod("append", signature(x = "whitneySquare", values = "whitneySquare"), function(x, values, after = length(x)) {
    return(whitneySquare(append(x@x, values@x, after), append(x@y, values@y, after), append(x@w, values@w, after)))
})

quadtree <- setClass("quadtree",
    slots = c(
        squares = "whitneySquare",
        # root = "numeric",
        children = "matrix",
        rep_x = "numeric",
        rep_y = "numeric"
    )
)

setMethod("initialize", "quadtree",
    function(.Object, squares, children, rep_x, rep_y) {
        .Object@squares <- squares
        .Object@children <- children
        .Object@rep_x <- rep_x
        .Object@rep_y <- rep_y
        return(.Object)
    }
)

# [i,j] = whether the i-th square contains the j-th point
belong_ <- function(squares, x, y) {
    outer(squares@x, x, `<=`) & outer(squares@x + squares@w, x, `>`) & outer(squares@y, y, `<=`) & outer(squares@y + squares@w, y, `>`)
}

# z shaped
quadsect_ <- function(square) {
    square@w <- c(c(0.5,0.5,0.5,0.5) %*% t(square@w))
    square@x <- c(c(1,1,1,1) %*% t(square@x)) + square@w * c(0,1,0,1)
    square@y <- c(c(1,1,1,1) %*% t(square@y)) + square@w * c(0,0,1,1)
    return(square)
}

partition_ <- function(x, y) {
    squares <- whitneySquare(0, 0, 1)
    children <- matrix(NA, 1, 4)
    rep_x <- c(NA_integer_)
    rep_y <- c(NA_integer_)

    q <- c(1)

    while (length(q) > 0) {
        queue <- squares[q]
        contain_matrix <- belong_(queue, x, y)
        contain_count <- rowSums(contain_matrix)

        contain_one <- contain_count == 1
        is_internal <- contain_count > 1

        sub_squares <- quadsect_(queue[is_internal])
        sub_contain_matrix <- belong_(sub_squares, x, y)
        replaced <- c()
        i <- 0
        j <- 0
        for (internal in is_internal) {
            j <- j + 1
            if (!internal) next

            i <- i + 1
            is_empty <- rowSums(sub_contain_matrix[(4*i-3):(4*i),]) == 0
            if (sum(is_empty) == 3) {
                squares[q[j]] <- sub_squares[4*i - 4 + which(!is_empty)]
                replaced <- append(replaced, q[j])
                is_internal[j] <- FALSE
            }
        }

        if (sum(contain_one) == 1) {
            rep_x[q[contain_one]] <- x[contain_matrix[contain_one,]]
            rep_y[q[contain_one]] <- y[contain_matrix[contain_one,]]
        } else {
            rep_index <- apply(contain_matrix[contain_one,],1,function(x) which(x)[1])
            rep_x[q[contain_one]] <- x[rep_index]
            rep_y[q[contain_one]] <- y[rep_index]
        }
        
        rep_x <- append(rep_x, rep(NA_integer_, 4*sum(is_internal)))
        rep_y <- append(rep_y, rep(NA_integer_, 4*sum(is_internal)))
        
        if (any(is_internal)) {
            i <- apply(matrix(contain_matrix[is_internal,], ncol=length(x)), 1, function(row) min(which(row)))
            rep_x[q[is_internal]] <- x[i]
            rep_y[q[is_internal]] <- y[i]
        }

        squares <- append(squares, quadsect_(queue[is_internal]))

        children[q[is_internal],] <- matrix(nrow(children) + seq_len(4*sum(is_internal)), ncol=4, byrow = TRUE)
        q <- c(replaced, nrow(children) + seq_len(4*sum(is_internal)))

        children <- rbind(children, matrix(NA, 4*sum(is_internal), 4))
    }

    return(quadtree(
        squares,
        children,
        rep_x,
        rep_y
    ))
}

# TODO vectorize this function, at least for u
ws_pairs <- function(tree, u, v, s) {
    is_leaf <- function(i) {
        return(all(is.na(tree@children[i,])))
    }

    is_well_separated <- function(i, j) {
        if (is_leaf(i) & is_leaf(j)) return(TRUE)
        distance <- max(abs(tree@rep_x[i] - tree@rep_x[j]), abs(tree@rep_y[i] - tree@rep_y[j]))
        radius <- tree@squares@w[i] + tree@squares@w[j]

        return(radius < s * distance)
    }

    if ((is_leaf(u) & is_leaf(v) & u == v)) return(NULL)
    if (is.na(tree@rep_x[u]) | is.na(tree@rep_x[v])) return(NULL)
    
    else if (is_well_separated(u, v)) {
        return(c(u, v))
    } else {
        if (tree@squares@w[u] < tree@squares@w[v]) {
            children <- tree@children[v,]
            v <- u
        } else {
            children <- tree@children[u,]
        }
        children <- Filter(Negate(is.null), children)
        return(unlist(sapply(children, function(child) ws_pairs(tree, child, v, s))))
    }
}

rect.whitney <- function(squares, x, y) {
    for (i in seq_along(squares)) {
        square <- squares[i]

        hasPoint <- any((square@x <= x & x < square@x + square@w) & (square@y <= y & y < square@y + square@w))
        col <- ifelse(hasPoint, "pink", "lightgreen")

        rect(square@x, square@y, square@x + square@w, square@y + square@w, col = col)
    }

    points(x, y, col = "red")
}

tree.debug <- function(tree) {
    is_leaf <- function(i) {
        return(all(is.na(tree@children[i,])))
    }

    for (i in seq_along(tree@squares)) {
        if (!is_leaf(i)) next

        square <- tree@squares[i]
        rep_x <- tree@rep_x[i]
        rep_y <- tree@rep_y[i]

        col <- ifelse(!is.na(rep_x), "pink", "lightgreen")

        rect(square@x, square@y, square@x + square@w, square@y + square@w, col = col)

        if (!is.na(rep_x)) {
            points(rep_x, rep_y, col = "red")
        }
    }

    data <- cbind(tree@squares@x, tree@squares@y, tree@squares@w, tree@children, tree@rep_x, tree@rep_y)
    colnames(data) <- c("x", "y", "w", "c1", "c2", "c3", "c4", "rep_x", "rep_y")
    View(data)
}

set.seed(2020)
### ANCHOR tree

n_arr <- seq(100,800,50)
mean_runtime <- c()
for (n in n_arr) {
    x <- runif(n)
    y <- runif(n)
    cat("Running with ", n, "\n")
    t <- as.numeric(summary(microbenchmark(partition_(x, y), times=10))["mean"])
    mean_runtime <- append(mean_runtime, t)
    cat("Result: ", t, "\n\n")
}
n_log_n <- (n_arr * log(n_arr))
plot(n_log_n, mean_runtime)
print(lm(formula = mean_runtime ~ n_log_n))
n <- 100:800
lines((n)*log(n), 0.08834*(n)*log(n)-52)

# n <- 800
# x <- runif(n)
# y <- runif(n)
# tree <- partition_(x, y)
# View(ws_pairs(tree, 1, 1, 0.01))
# plot(c(0,1), c(0,1), type = "n")
# tree.debug(tree)
# points(x, y, pch = 3)
# tree@squares[c(195,72)]
# rect.whitney(tree@squares[c(195,72)], x, y)
###

### ANCHOR plot
# range <- seq(100, 1000, 100)
# out <- range
# j <- 1
# for(i in range) {
#     x <- runif(i)
#     y <- runif(i)

#     tree <- partition_(x, y)

#     out[j] <- length(ws_pairs(tree, 1, 1, 0.5))/2
#     print(out[j])
#     j <- j + 1
# }

# plot(range, out)
###