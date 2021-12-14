devtools::load_all("../rrinterp", export_all = FALSE)

x <- c(4.723026, 5.367131, 5.733519, 6.401190, 6.492630, 6.509043, 6.567467, 7.159588, 7.701660, 8.284977, 8.810745)
y <- c(-5, 0.2259038, 4.4231622, 2.7695370, 3.9995587, 4.3818902, 3.5469489, 2.3960343, 1.1826827, 3.2680348, 5)
tau <- 5

interpolation <- rrinterpolate(x, y, -tau, tau)
rrplot(interpolation, x = x, y = y, limits = c(-tau, tau), autodiff = TRUE)
