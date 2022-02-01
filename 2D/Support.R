searchSupport.whitney <- function(squares, x, y) {
    if (length(x) != 1 || length(y) != 1) stop("searchSupport.whitney only supports singular point")
    base <- search.whitney(squares, x, y)

    if (is.na(base)) return(whitneySquare())
    if ((base@x + 0.2*base@w <= x & x < base@x + 0.8*base@w) & (base@x + 0.2*base@w <= y & y < base@y + 0.8*base@w)) return(base)

    # 10 9 8 7
    # 11     6
    # 12     5
    # 1  2 3 4
    x_key <- base@x + base@w/4 * c(-1, 1, 3, 5, 5, 5, 5, 3, 1, -1, -1, -1)
    y_key <- base@y + base@w/4 * c(-1, -1, -1, -1, 1, 3, 5, 5, 5, 5, 3, 1)

    squares <- search.whitney(squares, x_key, y_key, na.rm = TRUE)

    contain <- ((squares@x - 0.1*squares@w <= x & x < squares@x + 1.1*squares@w) & (squares@y - 0.1*squares@w <= y & y < squares@y + 1.1*squares@w))

    return(append(base, squares[contain]))
}

patch.whitney <- function(squares, x, y) {
    if (length(x) != 1 || length(y) != 1) stop("patch.whitney only supports singular point")
    p <- squares@field %*% c(1, x, y)

    delta_x <- (x - squares@x) / squares@w - 0.5
    delta_y <- (y - squares@y) / squares@w - 0.5

    phi_naught <- function(delta) {
        delta <- delta / 1.1
        return(-6*delta^5 + 15*delta^4 - 10*delta^3 +1)
    }

    phi <- phi_naught(delta_x) * phi_naught(delta_y)

    return(sum(p * phi / sum(phi)))
}