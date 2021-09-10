# Constrained-Interpolation
A 2021 summer research project

# Todo
- [ ] Check input.
- [x] Improve `as.piecewisePolynomial` for class `patching`. Currently, it results in redundant zero length segments, which could hinder the performance.
- [ ] Documentation.
- [ ] Turn it into a package.

# Known Issues
- [ ] The algorithm of `combinePiecewisePolynomial` is too naive in the sense that for a segment where only one piecewise polynomial is define, the result for that segment is just that polynomial. This doesn't make sense for operator such as `-` because we would expect the sign of right hand side to be flipped.