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
tau <- 3

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

threePointSolver.restricted <- function(data, tau. = tau) {
    slopes <- findSlope.beta.threePoints.restricted(data, tau = tau.)
    result <- interpolate.patch.onePointSlope(data, slopes, function(data, slope, tau.. = tau.) restrictedRange(data, slope, tau..))
    return(result)
}

my_plot <- function(interpolation, interval, data, usedual) {
    if (missing(usedual)) {
        f <- as.piecewisePolynomial(interpolation)
        if (is.null(f)) {
            usedual <- TRUE
        } else {
            interpolation <- f
            usedual <- FALSE
        }
    }

    if (usedual) {
        int <- dual(c(interval, rep(1, length(interval))), degree=2, length=length(interval), bydegree=TRUE)
        out <- predict(interpolation, int)
        df <- data.frame(
            x = interval,
            y = out[[,0]],
            y_prime = out[[,1]],
            y_prime2 = 2 * out[[,2]]
        )
    } else {
        d <- differentiate(interpolation)
        df <- data.frame(
            x = interval,
            y = predict(interpolation, interval),
            y_prime = predict(d, interval),
            y_prime2 = predict(differentiate(d), interval)
        )
    }
    df.melt <- melt(df, id = "x")
    p <- ggplot(df.melt, aes(x = x, y = value)) + 
        geom_line(aes(color = variable)) + 
        facet_grid(rows = variable ~ ., scales = "free_y")
    if (!missing(data)) {
        p <- p + geom_point(data = cbind(as.data.frame(data), variable="y"), 
             mapping = aes(x = x, y = y), 
             size = 1)
    }
    plot(p)
}

interpolation <- interpolate.patch.threePoint(data, threePointSolver.restricted)

int <- seq(min(x),max(x),0.05)
my_plot(interpolation, int, data)
