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

#' Minimize Quadratic with Linear Constraints
#' 
#' Minimizes `x^T A x + x^T B` under equality constraint `E x = F` and inequality constraint `G x >= H`
#' 
#' @details Uses Karush-Kuhn-Tucker method internally.
#' 
#' @param A A symmetric matrix describing the quadratic coefficients
#' @param B A vector describeing the linear coefficients
#' @param E A matrix for equality constraint
#' @param F A vector for equality constraint
#' @param G A matrix for inequality constraint
#' @param H A vector for inequality constraint
#' @param tol Tolerance. Number whose absolute value is less than `tol` is considered 0.
#' @return A list containing:
#' `x`: the desired `x`, see main description
#' `value`: the minimum value of `x^T A x + x^T B`
#' `lambda`: (optional) the Lagrangian multipliers for the equality constraints
#' `lambda`: (optional) the Lagrangian multipliers for the inequality constraints
#' @examples
#' A <- matrix(c(1, 0, 0, 2), nrow = 2, ncol = 2)
#' B <- c(-4, -4)
#' G <-matrix(c(-1, 1, -4, -1), nrow = 2, ncol = 2)
#' H <- c(-3, 0)
#' result <- min.quadratic(A, B, G = G, H = H)
min.quadratic <- function(A, B, E = NULL, F = NULL, G = NULL, H = NULL, tol = sqrt(.Machine$double.eps)) {
    if (ncol(A) < 1) stop("Dimension of `A` must be at least 1")
    if (ncol(A) != nrow(A)) stop("`A` must be a square matrix")
    if (ncol(A) != length(B)) stop("Dimensions of `A` and `B` don't match")

    if (is.null(E)) E <- matrix(nrow = 0, ncol = ncol(A))
    if (ncol(A) != ncol(E)) stop("Dimensions of `A` and `E` don't match")

    if (is.null(F)) F <- rep(0, nrow(E))
    if (nrow(E) != length(F)) stop("Dimensions of `E` and `F` don't match")

    if (is.null(G))  {
        M <- rbind(
            cbind(-2*A, t(E)),
            cbind(E, matrix(0,nrow=nrow(E),ncol=nrow(E)))
        )

        X <- limSolve::Solve(M, c(B, F), tol = tol)
        if (any(abs(M %*% X - c(B, F)) > tol)) return(NA)
        x <- X[seq_len(length(B))]
        value <- (t(x) %*% A %*% x + x %*% B)[1,1]

        result <- list(x=x, value=value)
        if (length(F) >= 1) result$lambda <- X[(length(B)+1):length(X)]
        return(result)
    }
    if (ncol(A) != ncol(G)) stop("Dimensions of `A` and `G` don't match")

    if (is.null(H)) H <- rep(0, nrow(G))
    if (nrow(G) != length(H)) stop("Dimensions of `G` and `H` don't match")

    result <- list(value=NA)

    activeInequality <- data.matrix(expand.grid(rep(list(0:1), nrow(G))))
    colnames(activeInequality) <- NULL

    for (i in seq_len(nrow(activeInequality))) {
        I <- which(activeInequality[i,] == 1)
        activeG <- G[I, ]
        activeH <- H[I]

        res <- min.quadratic(A, B, rbind(E, activeG), c(F, activeH))
        if (is.na(res)[1]) next

        mu <- if (length(I)==0) numeric() else res$lambda[(length(F)+1):length(res$lambda)]

        if (any(mu < 0)) next
        inactiveI <- which(activeInequality[i,] == 0)
        if (any(G[inactiveI,] %*% res$x < H[inactiveI])) next

        if (is.na(result$value) | res$value < result$value) {
            result <- res
            result$lambda <- if (length(F)==0) NULL else result$lambda[seq_len(length(F))]

            fullMu <- rep(0, nrow(G))
            fullMu[I] <- mu
            result$mu <- fullMu
        }
    }
    
    return(result)
}

#' Solving for beta
#' 
#' Find `beta` such that `t(beta) %*% A %*% beta + ||beta||_1` is minimzied subject to `B %*% beta = b`. `||beta||_1` is the 1-norm of `beta`
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
#' result <- solve.beta(A, B, b)
solve.beta <- function(A, B, b, tol = sqrt(.Machine$double.eps)) {
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

        result <- min.quadratic(Ahat, rep(1,6), Bhat, b, diag(6), rep(0,6), tol = tol)
        if (is.na(result)[1]) next

        if (is.na(min) | result$value < min) {
            min <- result$value
            beta[(I-1) %% 6 + 1] <- result$x * ifelse(I > 6, -1, 1)
        }
    }

    return(beta)
}