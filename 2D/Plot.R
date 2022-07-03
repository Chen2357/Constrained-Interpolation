library(plotly)

source("2D/WhitneyField.R")
source("2D/WhitneySquare.R")
source("2D/QuadTree.R")
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

W <- partition.quadtree(field)

# Plotting
x_res <- 50
y_res <- 50

value <- matrix(nrow = y_res, ncol = x_res)

for (i in seq_len(x_res)) {
    for (j in seq_len(y_res)) {
        checkpoint <- c(1.2 * (i-1) / (x_res-1) - 0.1, 1.2 * (j-1) / (y_res-1) - 0.1)

        support <- searchSupport.whitney(leaves.quadtree(W), checkpoint[1], checkpoint[2])
        value[j, i] <- patch.whitney(support, checkpoint[1], checkpoint[2])
    }
}

fig <- plot_ly(z = value)
fig <- fig %>% add_surface()

fig
