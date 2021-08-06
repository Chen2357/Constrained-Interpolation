library(methods)
library(limSolve)
source("Polynomial.R")
source("Data.R")
source("SlopeFinder.R")
source("PiecewisePolynomial.R")
source("Interpolater.R")
source("QuadraticProgramming.R")

### Old Main.R
# x <- c(1,2,3,4,5)
# y <- c(1,0,2,1,3)

# data <- pointData(x,y)

# # Example slope finding algorithms:
# # findSlope.quadratic
# # findSlope.constrainedQuadratic
# slope <- findSlope.quadratic(data)
# print(slope)

# # Example interpolation algorithms:
# # interpolate.patchQuadratic
# # interpolate.patchLinear
# # interpolate.joinLinear
# interpolation <- interpolate.patchQuadratic(data, slope)
# options(digits = 4)
# print(as.character(interpolation))
# print(as.data.frame(interpolation))

# plot(interpolation, xlab="x", ylab="y", ylim=range(-2,3.5))
# lines(differentiate(interpolation), col="blue")
# points(data, col="red")

### Three Points nonnegative beta Method
x <- c(0, 0.2, 1)
y <- c(4, 0, 1)

data <- pointData(x, y)
slopes <- findSlope.beta.threePoints(data)
print(slopes)

interpolation <- interpolate.patch.onePointSlope(data, slopes, quadratic.point.slope.extrema)
options(digits = 4)
print(as.character(interpolation))
print(as.data.frame(interpolation))

plot(interpolation, interval=seq(0,1,0.005), xlab="x", ylab="y")
lines(differentiate(interpolation), col="blue")
points(data, col="red")
