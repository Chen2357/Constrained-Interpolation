library(methods)
library(limSolve)
library(reshape2)
library(ggplot2)
source("Polynomial.R")
source("Data.R")
source("SlopeFinder.R")
source("PiecewisePolynomial.R")
source("Interpolater.R")
source("QuadraticProgramming.R")

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

# interpolation <- interpolate.patch.onePointSlope(data, slopes, quadratic.point.slope.extrema)
interpolation <- interpolate.patch.threePoint(data, threePointSolver)

int <- seq(leftMostBound(interpolation),rightMostBound(interpolation),0.005)
my_plot(interpolation, interval=int)

# plot(interpolation, interval=int, xlab="x", ylab="y")
# # plot(interpolation, interval=seq(0,10,0.001), xlab="x", ylab="y", ylim=range(-2,3.5))

# options(digits = 4)
# print(as.data.frame(interpolation))

# d <- differentiate(interpolation)
# lines(d * 0.2 + 1, interval=int, col="blue")
# points(data, col="red")

# lines(differentiate(d) * 0.1 + 1, interval=int, col="purple")

# print(predict(differentiate(interpolation), point.x(data)))

# ## Finding the point of maximum second derivative
# ex <- piecewisePolynomial.extrema(differentiate(differentiate(interpolation)))
# i <- which.max(abs(point.y(ex)))
# points(x=point.x(ex[i]), y=predict(interpolation, point.x(ex[i])), col = "purple")

# print(point.y(ex[i]))