# Constrained-Interpolation
A 2021 Summer to present research project

# Major Advancements
- [ ] Documentation.
- [x] Turn it into a package. ([info](https://swcarpentry.github.io/r-novice-inflammation/08-making-packages-R/))
    - See [rrinterp](https://github.com/Chen2357/rrinterp).

# Known Bugs
- [ ] The method `as.piecewisePolynomial` for class `patching` does not work as desired.
    - The issue is that polynomial coefficient blows up really quick. A potential solution is to add x-shift to polynomial.
- [ ] If the y value matches tau exactly, the program returns an error. `Error in L_inv %*% beta : non-conformable arguments`

# Code Improvements
- [x] Check user input.
- [x] Specify the number of nodes instead of the step size for the resolution.
- [x] Sort data by x value.
- [x] Implement operators where one side is missing.
- [x] Extract [Main.R](Main.R) to another file.
- [x] Remove `description` slot in `patchingFunction` because it is excessive.
- [ ] Implement `as.character` and printing for the `patching` class.
- [x] Check for zero length ranges in `setValidity` for `piecewisePolynomial`.

# Algorithm Improvements
- [ ] Check if simple quadratic can fit three points within restriction. If it does, use the quadratic curve.
    - This solves the issue where the derivative gets larger, even if the points are in a straight line.
- [ ] Experiment more with `solve.beta.cholesky`, which could be potentially faster than KKT, but currently, it does not behave like KKT.
- [ ] Use parallel computing to perform KKT. ([info](https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html))
- [ ] Add damping to `restrictedRange`.

# Minor Issues
- [ ] The algorithm of `combinePiecewisePolynomial` is too naive in the sense that for a segment where only one piecewise polynomial is define, the result for that segment is just that polynomial. This doesn't make sense for operator such as `-` because we would expect the sign of right hand side to be flipped.
- [ ] Dual number sometimes become `NaN` during calculation with the bump function.