quadtree <- setClass("quadtree",
    slots = c(
        squares = "whitneySquare",
        root = "numeric",
        children = "matrix"
    )
)

setMethod("initialize", "quadtree",
    function(.Object) {
        .Object@squares <- whitneySquare(0,0,1)
        .Object@root <- 1
        .Object@children <- matrix(NA_integer_,nrow = 1, ncol = 4)

        # validObject(.Object)
        return(.Object)
    }
)

bisect.quadtree <- function(tree, index) {
    bisection <- bisect.whitney(tree@squares[index])
    tree@children[index,] <- length(tree@squares) + seq_len(4*length(index))
    tree@squares <- append(tree@squares, bisection)
    tree@children <- rbind(tree@children, matrix(NA_integer_, nrow = 4*length(index), ncol = 4))
    return(tree)
}

partition.quadtree <- function(field) {
    tree <- quadtree()
    q <- c(1)

    x <- field@x
    y <- field@y

    while (length(q) > 0) {
        # This part need to be modified for readability
        queue <- tree@squares[q]
        contain_matrix <- outer(queue@x - queue@w, field@x, `<=`) & outer(queue@x + 2 * queue@w, field@x, `>`) & outer(queue@y - queue@w, field@y, `<=`) & outer(queue@y + 2 * queue@w, field@y, `>`)

        contain_indices <- apply(contain_matrix, 1, function(x) which(x)[1])

        contain_count <- rowSums(contain_matrix)
        contain_zero <- contain_count == 0
        contain_one <- contain_count == 1

        ok <- contain_one | contain_zero

        queue@field[contain_one,] <- field@coef[contain_indices[contain_one],]

        queue@field[!ok,] <- matrix(field@coef[contain_indices[!ok],], ncol = 3)

        tree@squares[q] <- queue
        len <- length(tree@squares)
        tree <- bisect.quadtree(tree, q[!ok])
        q <- len + seq_len(4*sum(!ok))
    }

    return(tree)
}

rect.quadtree <- function(tree, x, y) {
    for (i in seq_along(tree@squares)) {
        if (!all(is.na(tree@children[i,]))) next

        square <- tree@squares[i]

        hasPoint <- any((square@x <= x & x < square@x + square@w) & (square@y <= y & y < square@y + square@w))
        col <- ifelse(hasPoint, "pink", "lightgreen")

        rect(square@x, square@y, square@x + square@w, square@y + square@w, col = col)

        # TODO - omit this
        text(square@x + square@w/2, square@y + square@w/2, labels = square@field[1,1])
    }

    points(x, y, col = "red")
}

leaves.quadtree <- function(tree) {
    leaves <- apply(is.na(tree@children), 1, function(x) all(x))
    return(tree@squares[leaves])
}