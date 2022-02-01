# THIS FILE IS NOW OBSOLETE

# A comparison function is set on the whitneySquare.
# They are first compared by the `y` slot, and then the `x` slot.
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

# A whitneyDecomposition instance stores a whitneySquare instance along with a B-tree structure to allow fast searching.
# 
# `squares` is a whitneySquare instance that stores a collection of squares.
# `order` is the order of the B-tree. It should not be modified after initialization. (it is 5 by default)
# `root` is the index of the root note of the B-tree.
# `nodes` stores the information about each node of the B-tree. It is a matrix whose element nodes[x,y] is the y-th key (index of a square in `squares`) in the x-th node of the B-tree. Missing keys are represented by NA_integer_.
# `children` stores the information about the children of each node of the B-tree. It is a mtrix whose element children[x,y] is the y-th child of the x-th node of the B-tree. Missing children are represented by NA_integer_.
whitneyDecomposition <- setClass("whitneyDecomposition",
    slots = c(
        squares = "whitneySquare",
        order = "numeric",
        root = "numeric",
        nodes = "matrix",
        children = "matrix"
    )
)

setMethod("initialize", "whitneyDecomposition",
    function(.Object) {
        .Object@squares <- whitneySquare()
        .Object@order <- 5
        .Object@root <- NA_integer_
        .Object@nodes <- matrix(nrow = 0, ncol = .Object@order-1)
        .Object@children <- matrix(nrow = 0, ncol = .Object@order)

        # validObject(.Object)
        return(.Object)
    }
)

# `insert.whitney` is a B-tree insertion algorithm that allows the adding of a collection of square represented by a whitneySquare instance into a whitneyDecomposition instance.
insert.whitney <- function(decomposition, squares) {
    lastIndex <- length(decomposition@squares)
    decomposition@squares <- append(decomposition@squares, squares)

    order <- decomposition@order

    for (squareIndex in seq_along(squares)) {
        # `n` is the index of the node currently being searched.
        n <- decomposition@root
        # `key` is the index of the newly added square in `decomposition@squares`
        key <- lastIndex + squareIndex

        # `is.na(n)` is true when `decomposition` is empty with no squares.
        # in this case, we add the key of the square into a new node, and make it the root node.
        if (is.na(n)) {
            decomposition@nodes <- rbind(decomposition@nodes, c(key, rep(NA_integer_, order - 2)))
            decomposition@children <- rbind(decomposition@children, rep(NA_integer_, order))
            decomposition@root <- nrow(decomposition@nodes)
            next
        }

        # `nodePath` is a vector of the index of the nodes being searched.
        nodePath <- c()
        # `insertPath` is a vector of the index of the child node being searched.
        insertPath <- c()
        # `lastNotFull` is the last level being searched with missing element, the root node is on level 1.
        lastNotFull <- 0

        while (TRUE) {
            # `s` is a vector of the keys in node `n`
            s <- decomposition@nodes[n,]
            # `i` is the index of the child of node `n` that will be searched. (the index starts from 0)
            i <- sum(decomposition@squares[s[!is.na(s)]] < squares[squareIndex], na.rm=TRUE)
            # `notFull` is whether node `n` has missing element (which can potentially offer a place to store additional keys)
            notFull <- any(is.na(s))

            nodePath <- c(nodePath, n)
            insertPath <- c(insertPath, i)
            
            if (notFull) lastNotFull <- length(nodePath)
            
            # `is.na(decomposition@children[n,1])` is true when we are on the leaf node.
            if (is.na(decomposition@children[n,1])) {
                if (notFull) {
                    # if the leaf node is not full, we simply add the key to the missing slot.
                    decomposition@nodes[n,] <- append(s[1:(order-2)], key, i)
                    # `lastNotFull` is set to `NULL` to indicate that we don't need to do further action.
                    lastNotFull <- NULL
                }
                break
            } else {
                # if we are not on the left node, set the child node as the next node being searched.
                n <- decomposition@children[n,i+1]
            }
        }

        # `!is.null(lastNotFull)` is true when the leaf node is full, in which case, we need to take further actions.
        if (!is.null(lastNotFull)) {
            medium <- ceiling(order/2)
            # `hold` is the key that we are passing up the B-tree
            hold <- key
            # `holdChild` is the index of the right child after we split a node
            holdChild <- NA_integer_

            # `level` runs a backward loop from leaf level to the last level that is not full.
            for (level in length(nodePath):lastNotFull) {
                n <- nodePath[level]
                i <- insertPath[level]

                if (level != lastNotFull) {
                    # if we are not on the last level that is not full (meaning that node `n` is full), then we make an extended node (and its children) by inserting `hold` and `holdChild` into the insertion index we came from.
                    extendedNode <- append(decomposition@nodes[n,], hold, i)
                    extendedChildren <- append(decomposition@children[n,], holdChild, i+1)

                    # we split node `n` by replacing node `n` with the its first half, and create a new node from the second half. 
                    # we also hold on to the medium element that does not belong to either half, passing it up the B-tree.
                    decomposition@nodes[n,] <- c(extendedNode[1:(medium-1)], rep(NA_integer_, order-medium))
                    decomposition@children[n,] <- c(extendedChildren[1:medium], rep(NA_integer_, order-medium))

                    newNode <- extendedNode[(medium+1):order]
                    newNode <- c(newNode, rep(NA_integer_, medium - 1))
                    newChildren <- extendedChildren[(medium+1):(order+1)]
                    newChildren <- c(newChildren, rep(NA_integer_, medium - 1))

                    decomposition@nodes <- rbind(decomposition@nodes, newNode)
                    decomposition@children <- rbind(decomposition@children, newChildren)

                    hold <- extendedNode[medium]
                    holdChild <- nrow(decomposition@nodes)
                } else if (level == 0) {
                    # if we passed the root node (root node is on level 1), then we have to create a new node with what we are holding, and set it as the new root node.
                    newNode <- c(hold, rep(NA_integer_, order - 2))
                    newChildren <- c(decomposition@root, holdChild, rep(NA_integer_, order - 2))

                    decomposition@nodes <- rbind(decomposition@nodes, newNode)
                    decomposition@children <- rbind(decomposition@children, newChildren)

                    decomposition@root <- nrow(decomposition@nodes)
                } else {
                    # if we reached the last node that is not full, then we simply add what we are holding to the next missing slot in the node.
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

# `partition.whitney` is a B-tree insertion algorithm that allows the adding of a collection of square represented by a whitneySquare instance into a whitneyDecomposition instance.
partition._whitney <- function(field) {
    squares <- whitneySquare()
    queue <- whitneySquare(0,0,1)

    x <- field@x
    y <- field@y

    while (length(queue) > 0) {
        # This part need to be modified for readability
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

    return(insert.whitney(whitneyDecomposition(), squares))
}

# `rect.whitney` plots the a whitneyDecomposition instance, and indicates whether the each square contains a point.
rect._whitney <- function(decomposition, x, y) {
    for (i in seq_along(decomposition@squares)) {
        square <- decomposition@squares[i]

        hasPoint <- any((square@x <= x & x < square@x + square@w) & (square@y <= y & y < square@y + square@w))
        col <- ifelse(hasPoint, "pink", "lightgreen")

        rect(square@x, square@y, square@x + square@w, square@y + square@w, col = col)

        # TODO - omit this
        text(square@x + square@w/2, square@y + square@w/2, labels = square@field[1,1])
    }

    points(x, y, col = "red")
}

# FIXME - When a point is to the left of a square, it is not clear whether the point is in a square greater than less than
# `search.whitney` is a B-tree search algorithm that finds the squares that contain each of the points.
search._whitney <- function(decomposition, x, y, na.rm = FALSE) {
    result <- whitneySquare()

    squares <- W@squares
    for (i in seq_along(x)) {

        # FIXME - temporary O(n) algorithm
        contain <- (squares@x <= x[i] & x[i] < squares@x + squares@w) & (squares@y <= y[i] & y[i] < squares@y + squares@w)
        n <- which(contain)[1]
        if (is.na(n)) {
            if (!na.rm) result <- append(result, whitneySquare(NA_integer_, NA_integer_, NA_integer_))
        }
        else {
            result <- append(result, squares[n])
        }
        next

        n <- decomposition@root
        searching <- TRUE
        while (searching) {
            if (is.na(n)) {
                # the point is not in any of the squares
                if (!na.rm) result <- append(result, whitneySquare(NA_integer_, NA_integer_, NA_integer_))
                searching <- FALSE
                next
            }

            # the keys in the nodes
            nodes <- decomposition@nodes[n, ]
            # the number of keys in the node
            count <- sum(!is.na(nodes))
            # whether the point is found to be inside of a square in this node or whether the point is expected to be in one of the child node
            notFound <- TRUE

            for (j in seq_len(count)) {
                square <- decomposition@squares[nodes[j]]
                if (y[i] < square@y) {
                    # the point belongs to a square that is before the j-th square, so we search the j-th child
                    n <- decomposition@children[n, j]
                    notFound <- FALSE
                    break
                } else if (y[i] < square@y + square@w) {
                    if (x[i] < square@x) {
                        # the point belongs to a square that is before the j-th square, so we search the j-th child
                        n <- decomposition@children[n, j]
                        notFound <- FALSE
                        break
                    } else if (x[i] < square@x + square@w) {
                        # we have found that the point is inside this square.
                        result <- append(result, square)
                        notFound <- FALSE
                        searching <- FALSE
                        break
                    }
                }
            }

            if (notFound) {
                # the point belongs to a square that is greater than any of the squares in this node, so we search the rightmost child
                n <- decomposition@children[n, count+1]
            }
        }
    }

    return(result)
}