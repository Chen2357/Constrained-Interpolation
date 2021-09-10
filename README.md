# Constrained-Interpolation
A 2021 summer research project

# Todo
- [ ] Check user input.
- [x] Improve `as.piecewisePolynomial` for class `patching`. Currently, it results in redundant zero length segments, which could hinder the performance.
- [ ] Documentation.
- [ ] Turn it into a package.
- [ ] Experiment more with `solve.beta.cholesky`, which could be faster than KKT, but currently, it does not behave like KKT.
- [ ] Remove `description` slot in `patchingFunction` because it is excessive.
- [ ] Implement `as.character` and printing for the `patching` class.
- [ ] Implement operators where one side is missing.
- [ ] Extract [Main.R](Main.R) to another file.

# Known Issues
- [ ] The algorithm of `combinePiecewisePolynomial` is too naive in the sense that for a segment where only one piecewise polynomial is define, the result for that segment is just that polynomial. This doesn't make sense for operator such as `-` because we would expect the sign of right hand side to be flipped.
- [ ] Dual number sometimes become `NaN` during calculation with the bump function.