source("Polynomial.R")
source("PiecewisePolynomial.R")
source("DualNumber.R")
source("Patching.R")

tau <- 1
x0 <- 0
k <- 2
b <- 0

p <- polynomial(c(b-k*x0, k))

mu <- k^2 / (tau - abs(b))
delta <- (tau - abs(b)) / abs(k)

if (abs(k) < sqrt(.Machine$double.eps) | delta >= 1) {
    result <- patching(
        list(
            polynomial(0),
            p,
            polynomial(0)
        ),
        x0 + c(-1,0,1),
        patch.fifthDegree
    )
} else {
    q <- mu / 4 * polynomial(c(x0^2,-2*x0,1))
    result <- patching(
        list(
            polynomial(-tau),
            p+q,
            p,
            p-q,
            polynomial(tau)
        ), 
        x0 + delta * c(-2*sqrt(2),-1,0, 1,2*sqrt(2)),
        patch.fifthDegree
    )
}

int <- seq(-3,3,0.05)
y <- predict(result, int)

plot(int, y, type = "l")