# n=5 equal-first-four-moments witness

## Decisive finding

Two admissible mixing laws can agree exactly on `m1..m4` yet reverse the 50-year ranking of buffering (`R`) and doubled connectivity (`C`) for `n=5`, at the fixed parameters `p=1/10`, `q=3/100`, `e0=1/50`, `e1=9/10`, and `c=1/10`.

The common moment vector is

`(m0,m1,m2,m3,m4) = (1/1, 1/10, 3/100, 337/18000, 773993/51840000)`.

At the sharp lower fifth-moment endpoint `m5=1855440127/149299200000`, exact propagation gives `P_R-P_C < 0` (-0.00283857955333). At the sharp upper endpoint `m5=1990575127/149299200000`, it gives `P_R-P_C > 0` (0.0023575756861). The full exact rational differences are retained in `n5_equal_quadruples.json`.

For any labelled four-patch subset, a particular hit pattern with k hits has
probability `E[theta^k (1-theta)^(4-k)]`. Its degree is four, so these laws
agree on every complete four-patch environmental-hit distribution, and on
all lower-order marginals. Their five-way joint-hit probabilities differ.
This does not assert equality of occupancy histories.

## Evidence and assumptions

For strict-interior `m0..m4`, Schur complements give

`L5 = [m3,m4] [[m1,m2],[m2,m3]]^-1 [m3,m4]^T`

and

`U5 = m4 - [m2-m3,m3-m4] [[m0-m1,m1-m2],[m1-m2,m2-m3]]^-1 [m2-m3,m3-m4]^T`.

These are sharp because the odd-degree Hausdorff existence criterion on `[0,1]` is exactly `Hx=(m[i+j+1]) >= 0` and `H1-x=(m[i+j]-m[i+j+1]) >= 0` for `i,j=0..2`; see [Curto & Fialkow (1991), Theorem 4.1 and Remark 4.2, p. 622](https://www.math.uh.edu/~hjm/v017n4/0603CURTO.pdf). All principal minors of both matrices are nonnegative at both endpoints; one determinant is zero at each active endpoint. Thus both moment vectors are represented by probability laws on `[0,1]`.

The one-year kernels were interpolated exactly as degree-5 polynomials in the hazard `theta`. The rational witness was selected by a float64 `matrix_power(K,50)` screen, then independently propagated with exact rational matrix exponentiation. Every exact endpoint kernel is entrywise nonnegative, row-stochastic, and has absorbing state zero.

## Scope and residual risk

The parent audit also matched exact iterative 50-step propagation to the
matrix-power results, and reconstructed three-atom laws at both endpoints.
The reconstructed moments matched through degree five within 1.1e-17;
direct atomic-law floating-point propagation reproduced both signs.

This is an existence witness conditional on the specified iid annual
mixing-law model. It is not a full optimization over the n=5 moment class,
does not establish the prevalence or biological relevance of reversals,
and is not a theorem for arbitrary patch counts. Fable's independent audit
is the next review step.

Reproduce with
`uv run --python 3.12 --with numpy python ideation/b_screen/n5_equal_quadruples.py`.
