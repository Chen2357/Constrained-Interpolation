setMethod("symmetrize", "matrix", function(x) {
    x[lower.tri(x)] <- x[upper.tri(x)] <- (x[lower.tri(x)] + x[upper.tri(x)])/2
    return(x)
})

#' Pseudo Cholesky Decomposition
#' 
#' Perform Cholesky decomposition with remainder, that is find `L` and `R` such that `x = L %*% t(L) + R` where there are as many 0 as possible in `R`.
#' 
#' @details The algorithm is greedy, meaning the first nonzero component in `R` is as far down as possible in the diagonal.
#' 
#' @param x A symmetric matrix
#' @param tol Tolerance. Number whose absolute value is less than `tol` is considered 0.
#' @return A list containing:
#' `L`: the lower triangular matrix of the Cholseky decomposition
#' `R`: the remainder matrix 
cholesky <- function(x, tol = sqrt(.Machine$double.eps)) {
    L <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    colnames(L) <- colnames(x)
    j <- 1
    rownames(L) <- rownames(L)
    for (i in seq_len(ncol(L))) {
        if (x[i, i] > tol) {
            L[, j] <- x[, i] / sqrt(x[i, i])
            x <- x - L[, j] %*% t(L[, j])
            j <- j + 1
        }
    }
    return(list(L = L, R = x))
}

#' Solving for beta with Pseudo Cholesky Decomposition
#' 
#' Find `beta` such that `t(beta) %*% A %*% beta + ||beta||_1` is minimzied subject to `B %*% beta = b`. `||beta||_1` is the 1-norm of `beta`
#' 
#' @details Legacy function, use `\link[solve.beta]{solve.beta}` instead.
#' 
#' @param A A 6x6 symmetric postive-semidefinite matrix
#' @param B A 6x6 matrix
#' @param b A vector of 6 numbers
#' @param tol Tolerance. Number whose absolute value is less than `tol` is considered 0.
#' @return the desired `beta`, see main description
#' @examples
#' A <- diag(c(1 / 2, 1 / 3, 1 / 4, 1, 2, 3))
#' B <- diag(c(2, 3, 6, 0, 0, 0))
#' b <- c(1, 1, 1, 0, 0, 0)
#' result <- solve.beta.cholesky(A, B, b)
solve.beta.cholesky <- function(A, B, b, tol = sqrt(.Machine$double.eps)) {
    min <- NA
    beta <- NA

    betaPostive <- data.matrix(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1))
    colnames(betaPostive) <- NULL
    betaPostive <- cbind(betaPostive, !betaPostive)

    fullAhat <- kronecker(matrix(c(1, -1, -1, 1), nrow = 2, ncol = 2), A)

    for (i in 1:64) {
        I <- which(betaPostive[i, ] != 0)
        Ahat <- fullAhat[I, I]
        Bhat <- cbind(B, -B)[, I]

        Ahat <- rbind(cbind(Ahat, rep(0.5, 6)), t(c(rep(0.5, 6),0)))
        
        M <- cholesky(Ahat, tol = tol)
        if (any(abs(as.vector(M$R)[1:48]) > tol)) next

        trueA <- M$L[1:6, 1:6]
        trueB <- -M$L[7, 1:6]
        result <- limSolve::lsei(trueA, trueB, Bhat, b, diag(6), rep(0, 6), tol = tol, verbose = FALSE)
        if (result$IsError) next

        value <- result$solutionNorm + M$R[7,7]
        if (is.na(min) | value < min) {
            min <- value
            beta[(I-1) %% 6 + 1] <- result$x * ifelse(I > 6, -1, 1)
        }
    }

    return(beta)
}
