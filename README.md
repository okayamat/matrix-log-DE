# matrix-log-DE

MATLAB implementations for matrix logarithm algorithms
based on the double-exponential (DE) formula.

## Algorithms

- `OriginalAlg1.m` implements the DE formula-based algorithm of Tatsuoka et al.
  [1, 2].
- `ImprovedAlg2.m` implements the norm-free algorithm based on the DE formula
  proposed in the accompanying paper.
- `ImprovedAlg3.m` implements the error-bound-based algorithm, also based on
  the DE formula, proposed in the accompanying paper; it uses the error bound
  to select the number of quadrature nodes.

The improved algorithms are described in the accompanying paper. The routines
accept a matrix `A` and its dimension `m`; see each function signature for the
remaining inputs and returned values.

## Numerical experiments

The following scripts run experiments and write tab-delimited `.dat` files:

- `PointError.m` compares the relative errors of the three algorithms as the
  number of quadrature nodes varies.
- `TimeError.m` records execution times and relative errors for the three
  algorithms as the maximum number of nodes varies.
- `ErrorBound.m` records the error bound returned by `ImprovedAlg3.m` for
  several maximum node counts.

The experiments use test matrices with condition numbers
`10`, `10^4`, and `10^7`. They construct each test matrix as
`A = Q*D*Q'`, where `Q` is loaded from the corresponding size directory
(e.g., `50/Qdata.csv`) and `D` is a positive diagonal matrix. The scripts
process matrix dimensions `m = 50, 100, 200, 400, 800, 1600` (the
`ErrorBound.m` script uses `m = 50` only).

Run the scripts in MATLAB with the current folder set to this directory.
The required `Qdata.csv` files must be available in the corresponding
numbered folders (`50`, `100`, `200`, `400`, `800`, and `1600`). These folders
contain data for matrices of the corresponding sizes and are also used for
output.
Ensure these folders exist before running the scripts.

Each output file contains two columns. `PointError.m` writes the number of
quadrature nodes and the relative error in the matrix 2-norm. `TimeError.m`
writes the elapsed time in seconds and the relative error in the matrix 2-norm.
`ErrorBound.m` writes the number of quadrature nodes and the error bound
`En` supplied by `ImprovedAlg3.m`. Output filenames identify the algorithm,
experiment, and condition number.

## References

[1] F. Tatsuoka, T. Sogabe, Y. Miyatake, and S.-L. Zhang, “Algorithms for
the Computation of the Matrix Logarithm Based on the Double Exponential
Formula,” *Journal of Computational and Applied Mathematics*, vol. 373,
article 112396, 2020.

[2] F. Tatsuoka, T. Sogabe, T. Kemmochi, and S.-L. Zhang, “On the
Convergence Rate of the Double-Exponential Formula for the Computation of
the Matrix Logarithm,” *RIMS Kōkyūroku*, vol. 2167, pp. 1--9,
2020 (in Japanese).
