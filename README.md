# Constrained-Interpolation
A 2021 summer research project

# Major Advancements
- [ ] Documentation.
- [ ] Turn it into a package. ([info](https://swcarpentry.github.io/r-novice-inflammation/08-making-packages-R/))

# Known Bugs
- [ ] The method `as.piecewisePolynomial` for class `patching` does not work as desired.
- [ ] If the y value matches tau exactly, the program returns an error. `Error in L_inv %*% beta : non-conformable arguments`

# Code Improvements
- [ ] Check user input.
- [ ] Specify the number of nodes instead of the step size for the resolution.
- [ ] Sort data by x value.
- [ ] Use parallel computing to perform KKT. ([info](https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html))
- [ ] Implement operators where one side is missing.
- [ ] Extract [Main.R](Main.R) to another file.
- [ ] Remove `description` slot in `patchingFunction` because it is excessive.
- [ ] Implement `as.character` and printing for the `patching` class.
- [ ] Check for zero length ranges in `setValidity` for `piecewisePolynomial`.

# Algorithm Improvements
- [ ] Check if simple quadratic can fit three points within restriction. If it does, use the quadratic curve.
    - This solves the issue where the derivative gets larger, even if the points are in a straight line.
- [ ] Experiment more with `solve.beta.cholesky`, which could be potentially faster than KKT, but currently, it does not behave like KKT.

# Minor Issues
- [ ] The algorithm of `combinePiecewisePolynomial` is too naive in the sense that for a segment where only one piecewise polynomial is define, the result for that segment is just that polynomial. This doesn't make sense for operator such as `-` because we would expect the sign of right hand side to be flipped.
- [ ] Dual number sometimes become `NaN` during calculation with the bump function.