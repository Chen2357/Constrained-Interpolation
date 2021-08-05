library(limSolve)
source("QuadraticProgramming.R")

A <- diag(c(1 / 2, 1 / 3, 1 / 4, 1, 2, 3))
B <- diag(c(2, 3, 6, 0, 0, 0))
b <- c(1, 1, 1, 0, 0, 0)
result <- solve.beta(A, B, b)
