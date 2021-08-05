library(limSolve)

#' Minimize Quadratic with Linear Constraints
#' 
#' Minimizes `x^T A x + x^T B` under equality constraint `E x = F` and inequality constraint `G x >= F`
#' 
#' @param A A symmetric matrix describing the quadratic coefficients
#' @param B A vector describeing the linear coefficients
#' @param E A matrix for equality constraint
#' @param F A vector for equality constraint
#' @param G A matrix for inequality constraint
#' @param H A vector for inequality constraint
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
min.quadratic <- function(A, B, E = NULL, F = NULL, G = NULL, H = NULL) {
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

        X <- solve(M, c(B, F))
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

        res <- try(min.quadratic(A, B, rbind(E, activeG), c(F, activeH)))

        mu <- if (i==0) numeric() else res$lambda[(length(F)+1):length(res$lambda)]

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
