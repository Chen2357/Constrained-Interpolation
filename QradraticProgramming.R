library(limSolve)

#' Pseudo Cholesky Decomposition
#' 
#' Perform Cholesky decomposition with remainder, that is find `L` and `R` such that `x = L %*% t(L) + R`. 
#' 
#' @param x A symmetric matrix
#' @param tol Tolerance
#' @return A list containing:
#' `L`: the lower triangular matrix of the Cholseky decomposition
#' `R`: the remainder matrix 
cholesky <- function(x, tol = .Machine$double.eps^0.5) {
    L <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    colnames(L) <- colnames(x)
    rownames(L) <- rownames(L)
    for (i in seq_len(ncol(L))) {
        if (x[i, i] > tol) {
            L[, i] <- x[, i] / sqrt(x[i, i])
            x <- x - L[, i] %*% t(L[, i])
        }
    }
    return(list(L = L, R = x))
}

solve.beta <- function(A, B, b, tol = .Machine$double.eps^0.5) {
    min <- NA
    sol <- NA

    betaPlus <- expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
    betaMinus <- rep(1, 6) - expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
    colnames(betaPlus) <- NULL
    colnames(betaMinus) <- NULL

    betaNonzero <- cbind(betaPlus, betaMinus)
    colnames(betaNonzero) <- NULL
    betaNonzero <- data.matrix(betaNonzero)

    fullAhat <- kronecker(matrix(c(1, -1, -1, 1), nrow = 2, ncol = 2), A)

    for (i in 1:64) {
        I <- which(betaNonzero[i, ] != 0)
        Ahat <- fullAhat[I, I]
        Bhat <- cbind(B, -B)[, I]

        Ahat[lower.tri(Ahat)] <- Ahat[upper.tri(Ahat)] <-sym <- (Ahat[lower.tri(Ahat)] + Ahat[upper.tri(Ahat)])/2
        Ahat <- rbind(cbind(Ahat, rep(0.5, 6)), t(c(rep(0.5, 6),0)))
        
        M <- cholesky(Ahat, tol = tol)
        if (any(abs(as.vector(M$R)[1:48]) > tol)) next

        trueA <- M$L[1:6, 1:6]
        trueB <- -M$L[7, 1:6]
        result <- lsei(trueA, trueB, Bhat, b, diag(6), rep(0, 6), verbose = FALSE)
        if (result$IsError) next

        value <- result$solutionNorm + M$R[7,7]
        if (is.na(min) | value < min) {
            min <- value
            sol <- result$X
            I <- which(betaNonzero[i, 1:6] == 0)
            if (length(I)>0) sol[I] <- -sol[I]
        }
    }

    return(list(solution = sol, value = min))
}

A <- diag(c(1 / 2, 1 / 3, 1 / 4, 1, 2, 3))
B <- diag(c(2, 3, 6, 0, 0, 0))
b <- c(1, 1, 1, 0, 0, 0)
result <- solve.beta(A, B, b)