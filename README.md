# matrix-log-DE

MATLAB implementations and numerical data for matrix logarithm algorithms
based on the double-exponential (DE) formula.

## Algorithms

- `OriginalAlg1.m` implements the DE formula-based matrix logarithm
  algorithm developed by Tatsuoka et al. [1, 2].
- `ImprovedAlg2.m` implements a norm-free adaptive algorithm for computing the
  matrix logarithm with the DE formula.
- `ImprovedAlg3.m` implements an error-bound-based algorithm for the same
  problem; it uses the bound to select the number of quadrature nodes.

The latter two algorithms are described and analyzed in th accompanying
paper [3]. The routines accept a matrix `A` and its dimension `m`; see each
function signature for the remaining inputs and returned values.

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
Each matrix-size directory (`50`, `100`, `200`, `400`, `800`, and `1600`)
contains its corresponding `Qdata.csv` input file and is also used for output.
`ErrorBound.m` uses only the `50` directory.

Each output file contains two columns. `PointError.m` writes the number of
quadrature nodes and the relative error in the matrix 2-norm. `TimeError.m`
writes the elapsed time in seconds and the relative error in the matrix 2-norm.
`ErrorBound.m` writes the number of quadrature nodes and the error bound
`En` supplied by `ImprovedAlg3.m`. Output filenames identify the algorithm,
experiment, and condition number.

## Computational environment

The numerical experiments were performed on a MacBook Air with an Apple M4
chip and 24 GB of memory, running macOS Tahoe 26.6.2 and MATLAB R2024b.
The maximum number of computational threads reported by MATLAB's
`maxNumCompThreads` function was 10; this value was queried rather than
explicitly set. Execution times were measured using
MATLAB's `timeit` function.

## References

[1] F. Tatsuoka, T. Sogabe, Y. Miyatake, and S.-L. Zhang, "Algorithms for
the computation of the matrix logarithm based on the double exponential
formula," *Journal of Computational and Applied Mathematics*, vol. 373,
article 112396, 2020.

[2] F. Tatsuoka, T. Sogabe, T. Kemmochi, and S.-L. Zhang, "On the
convergence rate of the double-exponential formula for the computation of
the matrix logarithm," *RIMS Kôkyûroku*, vol. 2167, pp. 1--9,
2020 (in Japanese).

[3] T. Okayama and T. Setogawa, "Improved algorithms for computing the matrix
logarithm based on the double-exponential formula," arXiv:2610.25361, 2026.
