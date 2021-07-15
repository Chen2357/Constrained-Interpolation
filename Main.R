library(methods)
source("Polynomial.R")
source("Data.R")
source("SlopeFinder.R")
source("PiecewisePolynomial.R")
source("Interpolater.R")

x <- c(1,2,3,4,5)
y <- c(1,3,3,1,3)

data <- pointData(x=x,y=y)

# Example slope finding algorithms:
# quadraticSlopeFind
# constrainedQuadraticSlopeFind
slopeData <- constrainedQuadraticSlopeFind(data)
print(slope(slopeData))

# Example interpolation algorithms:
# quadraticPatchInterpolate
# linearPatchInterpolate
# quadraticJointInterpolate
interpolation <- quadraticPatchInterpolate(slopeData)
print(interpolation, format = "table", digits = 5)

plot(interpolation, xlab="x", ylab="y")
plot(data, col="red")
