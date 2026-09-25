# Independent review package: BWBrook/deep

Reviewed commit: `6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3`.

Start with **review.md**. This package contains independent reviewer code and
outputs, not a clone or modification of the private repository. No credentials
or signed download links are included. All source links in the report are pinned.

## Contents

- `review.md`: scientific assessment, ecological relevance, priority actions,
  verification scope, and literature-access qualifications.
- `independent_checks.py` / `.json`: fresh labelled transition model, exact
  propagation, reference laws, attained ranges, trajectories and shape witnesses.
- `independent_certificates.py` / `independent_certificates_vertical.json`:
  exact exchangeable-polygon global certificates at all 16 reference cases.
- `additional_checks.py` / `.json`: all 32 scalar curved boundaries,
  reference maximin-persistence certificates, equal-triples witnesses,
  unimodal cutoff and strictly positive-support examples.
- `validate_package.py` / `validation.json`: exact joins between certificates
  and attained values; comparison against all 64 rounded table endpoints.
- `*.log`: execution records. `environment.json` and `manifest.sha256` record
  the tested environment and package hashes.

## Reproduce

Tested here with CPython 3.13.5 and NumPy 2.3.5. The code uses Python integers
and `fractions.Fraction` for every claimed exact sign and bound. NumPy arrays
have object dtype and do not substitute floating-point arithmetic in proofs.
Decimal values in JSON are display conveniences; the `exact` strings govern.

Create a virtual environment, install the recorded dependency, and run:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python independent_checks.py
OPENBLAS_NUM_THREADS=1 python independent_certificates.py --vertical
OPENBLAS_NUM_THREADS=1 python additional_checks.py
python validate_package.py
```

Do not run with `python -O`: assertions are part of validation. The scripts
write receipts beside themselves. Rerunning changes recorded runtime fields;
therefore output file hashes may differ even when the exact scientific payload
is unchanged. The last script verifies the exact relationships independently
of those runtime metadata.

## Why these are class certificates

The finite-exchangeable four-patch feasible set is a polygon in third- and
fourth-order co-hit probabilities. A strictly positive vertical derivative of
Delta on a containing rectangle excludes interior extrema. Monotone polygon
edges then place the extrema at vertices. The rectangle can include invalid
probability laws: the polynomial extension remains defined there, and bounding
its derivative on a superset is legitimate.

For the scalar class, the feasible set is a curved lens within that polygon.
The same vertical derivative check moves extrema to its two curved boundaries.
Their degree-100 Bernstein coefficients lie between attained endpoint values.
This certifies endpoint extrema without assuming that every nonlinear objective
on a convex set is vertex-optimised. Polynomial evaluations at interior points
are checked against independently factorised transition propagation.

The reference maximin analysis separately certifies P_R and P_C; it does not
infer their extrema from Delta's extrema. It certifies nonzero vertical signs
and monotone edges for each objective at q=.03 and .06. The full-support
examples, in contrast, are just exact feasible witnesses, not class optimisations.

## Scope

This code primarily verifies the four-patch results. It is not a complete
reproduction of all repository analyses, a production general-purpose solver,
an empirical model, or a formal proof assistant. The review states what was
and was not independently executed. The larger-n and temporal-asymmetry results
are not certified by these scripts. A class label on an attained range in
`independent_checks.json` becomes a global claim only after joining it to the
appropriate certificate in the subsequent scripts.
