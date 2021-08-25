library(methods)
library(limSolve)
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

interpolation <- interpolate.patch.threePoint(data, threePointSolver)
plot(interpolation, xlab="x", ylab="y")
# plot(interpolation, interval=seq(0,10,0.001), xlab="x", ylab="y", ylim=range(-2,3.5))

options(digits = 4)
print(as.data.frame(interpolation))

lines(differentiate(interpolation), col="blue")
points(data, col="red")

## Finding the point of maximum second derivative
ex <- piecewisePolynomial.extrema(differentiate(differentiate(interpolation)))
i <- which.max(abs(point.y(ex)))
points(x=point.x(ex[i]), y=predict(interpolation, point.x(ex[i])), col = "purple")

print(point.y(ex[i]))