source("WhitneySquare.R")
source("WhitneyDecomposition.R")

x <- c(0.2,0.2,0.4,0.3)
y <- c(0.3,0.1,0.5,0.9)

W <- partition.whitney(x, y)

plot(c(0,1), c(0,1), type = "n")
rect.whitney(W, x, y)
