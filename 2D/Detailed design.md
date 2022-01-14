# 2D Interpolation

## Detailed design

### Whitney field
- `whitneyField` 
    - A collection of linear polynomials whose coefficients are stored in a $n\times 3$ matrix
    - `whitneySquare` will be extended to have reference to the corresponding polynomial in the `whitneyField`. We will need to assign each Whitney square with a polynomial. We can consider several implementations:
        - Storing the index (row number) of the corresponding polynomial in `whitneyField` (possibly also extend `whitneySquare` to store the `whitneyField`).
        - Directly store the linear polynomials in a matrix in `whitneyField`. In which case, we don't need the `whitneyField` class at all.

### Record the support of functions
- `supportRect`
    - An array of rectangle (with info about bottom right corner coordinates, height, and width).
    - Each rectangle is associated with a linear polynomial in `whitneyField` represented by the index (row number) of the polynomial in the matrix.
    - The initialization of `supportRect` will use `whitneySquare` where each square has support on its $1.1$ times dilation region.

- `supportMap`
    - A B-tree structure for `supportRect` similar for how `whitneyDecomposition` stores a B-tree structure for `whitneySquare`.
    - It can take a point and returns the indices of the polynomials (corresponding to each `whitneySquare`) that has support on the point.

### Whitney interpolation
- `whitneyInterpolation`
    - Contains `supportMap`, `whitneySquare` (`whitneyField`), patching method (perhaps represented by characters `"5"` or `"bump"` to signify which patching method to use).
    - When user feed in a point, it perform search in `supportMap` to get the `whitneySquare`s that have support on the point. Polynomial values are weighted by the patching function values, and the  interpolated value for the point will be returned.

### Plotting
- To be completed.

## Future directions

- Hyper-dual numbers for multivariate automatic differentiation.
