source("2D/WhitneyField.R")
source("2D/WhitneySquare.R")
# source("2D/WhitneyDecomposition.R")
source("2D/Support.R")

x <- c(0.2,0.2,0.4,0.3)
y <- c(0.3,0.1,0.5,0.9)
coef <- matrix(c(
    1,0,0,
    2,0,0,
    3,0,0,
    4,0,0
), byrow = TRUE, ncol = 3)

field <- whitneyField(x, y, coef)

W <- partition.whitney(field)

plot(c(0,1), c(0,1), type = "n")
rect.whitney(W, x, y)

# checkpoint <- c(0.37, 0.45)
# print(search.whitney(W, checkpoint[1], checkpoint[2]))

# squares <- W@squares
# contain <- (squares@x <= checkpoint[1] & checkpoint[1] < squares@x + squares@w) & (squares@y <= checkpoint[2] & checkpoint[2] < squares@y + squares@w)
# i <- which(contain)
# print(i)
# View(cbind(W@nodes, W@children))

checkpoint <- c(0.375, 0.4)
support <- searchSupport.whitney(W, checkpoint[1], checkpoint[2])
print(support)
print(patch.whitney(support, checkpoint[1], checkpoint[2]))