library(methods)
library(limSolve)
library(reshape2)
library(ggplot2)
source("Utilities.R")
source("Polynomial.R")
source("Data.R")
source("SlopeFinder.R")
source("PiecewisePolynomial.R")
source("Patching.R")
source("Interpolater.R")
source("QuadraticProgramming.R")
source("DualNumber.R")

x <- c(1,2,3,4,5)
y <- c(1,0,2,1,3)

data <- pointData(x,y)

threePointSolver <- function(data) {
    result <- quadraticPolynomial(data)
    ce <- coef(result)
    if (ce[2] * ce[2] - 4 * ce[1] * ce[3] > 0) {
        slopes <- findSlope.beta.threePoints(data)
        result <- interpolate.patch.onePointSlope(data, slopes, quadratic.point.slope.extrema)
    }
    return(result)
}

my_plot <- function(interpolation, interval = seq(leftMostBound(interpolation),rightMostBound(interpolation),0.05)) {
    d <- differentiate(interpolation)
    df <- data.frame(
        x = interval,
        y = predict(interpolation, interval),
        y_prime = predict(d, interval),
        y_prime2 = predict(differentiate(d), interval)
    )
    df.melt <- melt(df, id = "x")
    plot(
        ggplot(df.melt, aes(x = x, y = value)) + 
        geom_line(aes(color = variable)) + 
        facet_grid(rows = variable ~ ., scales = "free_y")
    )
}

interpolation <- interpolate.patch.threePoint(data, threePointSolver)
interpolation <- as.piecewisePolynomial(interpolation, min(x), max(x))

int <- seq(min(x),max(x),0.05)
my_plot(interpolation, interval=int)
