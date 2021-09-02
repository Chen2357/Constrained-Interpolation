# Constrained-Interpolation
A 2021 summer research project

# Todo
- [ ] Check input.
- [ ] Make patching using the bump function possible.
- [ ] Add points to the final plot.
- [ ] Implement [range.pdf](Resources/range/pdf).

# Known Issues
- [ ] The algorithm of `combinePiecewisePolynomial` is too naive in the sense that for a segment where only one piecewise polynomial is define, the result for that segment is just that polynomial. This doesn't make sense for operator such as `-` because we would expect the sign of right hand side to be flipped.