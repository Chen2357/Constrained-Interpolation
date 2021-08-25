# Constrained-Interpolation
A 2021 summer research project

# Todo
- [x] Implement [kkt.pdf](Resources/kkt.pdf)
- [x] Successfully interpolate three points with the beta method
- [ ] Implement interpolating n points.
- [x] Find the extrema points of polynomial and piecewise polynomial.
  - [ ] ~~Require the ranges in piecewise polynomial to be ordered.~~
- [x] Improve interface
  - [x] Documentation
  - [ ] ~~Possibly turn the project into a library with [devtools](https://www.rdocumentation.org/packages/devtools/versions/2.4.2) package~~ (dropped because there is no perceivable benefit)

# Known Issues
- [ ] The algorithm of `combinePiecewisePolynomial` is too naive in the sense that for a segment where only one piecewise polynomial is define, the result for that segment is just that polynomial. This doesn't make sense for operator such as `-` because we would expect the sign of right hand side to be flipped.