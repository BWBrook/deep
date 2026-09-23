# A local nonexchangeable witness with identical pairwise hit laws

Astra, 2026-09-23. This addresses the construction objection in correspondence
#470. It proves a limited existence result: exact exchangeability is not
necessary for the ranking reversal. It does not validate a spatial model,
calibrate heterogeneity, or extend the unimodal-class theorem.

## Coupling argument

Let Qminus and Qplus be the one-step distributions of the complete labelled
hit vector under the existing four-patch endpoint hazard laws. Their complete
pairwise marginals agree, and their action differences satisfy

    Delta_minus = -0.12052030700439428...,
    Delta_plus  = +0.03209273498500963... .

Keep the original mortality/recolonisation rules and the same buffered patch
under R. For any common hit-vector distribution H, set

    Q_s,epsilon = (1-epsilon) Q_s + epsilon H,  s in {minus,plus}.

The pairwise marginals of these two new distributions remain equal, since
marginalisation is linear. Draw independently from this fixed distribution
each model step, as in the original temporal specification.

Couple an original and a perturbed process by retaining the original shock
vector with probability 1-epsilon each step and using H otherwise. Couple
all mortality and recolonisation randomness as well. With probability
(1-epsilon)^T, no replacement occurs and both occupancy trajectories agree.
Therefore each action's persistence changes by at most
eta=1-(1-epsilon)^T. Subtraction gives

    |Delta_s,epsilon - Delta_s| <= 2 eta <= 2T epsilon.

If the original pair regret r=min(-Delta_minus,Delta_plus) exceeds 2 eta,
both opposite signs persist and the new pair regret is at least r-2 eta.
This argument uses labelled trajectories, so it does not require a valid
occupancy-count aggregation for the perturbed model.

## Explicit asymmetry, preserving the original p and q

For x in {0,1}^4 let k=sum(x), and write the Glo hit-vector law as

    Qminus(x) = (2/3) 1{k=0} + (1/3)(3/10)^k(7/10)^(4-k).

Choose

    H(x) = Qminus(x) + (-1)^(x1+x2+x3)/1600.

The signed term sums to zero. In any marginal on at most two patches, at
least one of x1,x2,x3 is summed out, so its contribution also cancels. Thus
H retains every original pairwise hit table: (00,01,10,11) has probabilities
(83,7,7,3)/100, giving p=1/10 and q=3/100.

H is a probability law: its minimum pattern probability is 83/40000>0.
It is nonexchangeable because H(1000)=1347/40000 differs from
H(0001)=1397/40000. Both perturbed laws inherit this asymmetry for every
epsilon>0, while retaining exactly the original common pairwise tables.
In particular, common marginal rates and common pairwise correlations do
not imply exchangeability of the full shock distribution.

Take epsilon=1/4000 and T=50. The original exact endpoint results give
r>3209/100000, while 2T epsilon=1/40. Therefore the perturbed pair has

    regret > 709/100000 = 0.00709,

or more than 0.709 percentage points of persistence, with opposite rankings.
Using 2 eta instead of its linear bound gives the stronger lower bound
0.007245249280097929... . Both are conservative bounds, not exact evaluations
of the perturbed action differences.

## Verification and interpretation

`python3 ideation/b_screen/asymmetric_pairwise_witness.py` checks all 120
pairwise-pattern equalities across the five laws, probability feasibility,
nonexchangeability of H and the perturbed pair, and the rational regret
inequalities. Its JSON retains the complete labelled laws, exact coupling
bound, original exact endpoint values and source hashes. The existing n4
endpoint certificate is reused, not replayed. No exchangeable transition
kernel is applied to an asymmetric shock law.

This is a small, artificial perturbation designed to test whether exact
symmetry is essential. It rules out that narrow objection. It does not show
that substantial spatial heterogeneity preserves the reported loss magnitude,
that these perturbations represent any landscape, or that the original
small-n full-class bounds cover nonexchangeable laws. Those questions require
a separate modelling argument.
