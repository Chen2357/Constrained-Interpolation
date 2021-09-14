rrplot <- function(interpolation, interval, x.point, y.point, limits, res = 1000, autodiff) {
    if (missing(interval)) {
        if (!missing(x.point)) {
            interval <- seq(min(x.point), max(x.point), (max(x.point)-min(x.point))/res)
        } else {
            stop("Unable to infer the interval to plot")
        }
    }

    if (!missing(x.point)) {
        if(!missing(y.point)) {
            data <- pointData(x, y)
        } else {
            data <- pointData(x, predict(interpolation, x))
        }
    }

    if (missing(autodiff)) {
        f <- as.piecewisePolynomial(interpolation)
        if (is.null(f)) {
            usedual <- TRUE
        } else {
            interpolation <- f
            usedual <- FALSE
        }
    }

    if (autodiff) {
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
    if (!missing(limits)) {
        p <- p + geom_hline(
            data = cbind(data.frame(limits), variable="y"),
            aes(yintercept = limits),
            color = "blue",
            linetype = "dashed"
        )
    }
    plot(p)
}

rrinterpolate.plot <- function(x, y, min, max, res = 1000, usedual) {
    interpolation <- rrinterpolate(x, y, min, max)
    interval <- seq(min(x), max(x), (max(x) - min(x)) / res)
    rrplot(interpolation, interval, pointData(x, y), usedual)
}