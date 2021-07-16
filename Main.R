library(methods)
source("Polynomial.R")
source("Data.R")
source("SlopeFinder.R")
source("PiecewisePolynomial.R")
source("Interpolater.R")

x <- c(1,2,3,4,5)
y <- c(1,0,2,1,3)

data <- pointData(x,y)

# Example slope finding algorithms:
# quadraticSlopeFind
# constrainedQuadraticSlopeFind
slope <- quadraticSlopeFind(data)
print(slope)

# Example interpolation algorithms:
# quadraticPatchInterpolate
# linearPatchInterpolate
# quadraticJointInterpolate
interpolation <- quadraticPatchInterpolate(data, slope)
options(digits = 4)
print(as.character(interpolation))
print(as.data.frame(interpolation))

plot(interpolation, xlab="x", ylab="y", ylim=range(-2,3.5))
lines(differentiate(interpolation), col="blue")
points(data, col="red")
