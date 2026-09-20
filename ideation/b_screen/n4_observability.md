# What a complete hit-count record can distinguish

This is a response to correspondence thread 49, message 401. The inverse of
a difference between rare-event probabilities is not a sample-size or
discrimination calculation. It cannot establish practical unobservability.

For the equal-triples witness in `n4_audit.md`, the two laws have the same
moments through degree three and different fourth moments. The complete
yearly hit count K has distribution

`P(K=k) = choose(4,k) E[theta^k (1-theta)^(4-k)]`.

Changing m4 at fixed lower moments changes these five probabilities in the
direction `(1,-4,6,-4,1)`. Thus a full record informs m4 through all hit
counts, not only through years with K=4. Exchangeability means patch labels
provide no extra information about which of these two laws generated a
given count. This statement concerns environmental hits, not occupancies.

| K | Lower-boundary law | Upper-boundary law |
|---|---:|---:|
| 0 | 0.713100 | 0.715843 |
| 1 | 0.223600 | 0.212629 |
| 2 | 0.030600 | 0.047057 |
| 3 | 0.015600 | 0.004629 |
| 4 | 0.017100 | 0.019843 |

## A deliberately specified finite-sample benchmark

Assume the two candidate laws are known, have equal prior probabilities,
and generate independent years. All four environmental hits are observed
without error. With zero-one loss, the minimum average probability of
selecting the wrong law is half the sum of the smaller of its two
likelihoods over all possible records. Since the two laws favour opposite
actions, this is also the minimum wrong-action probability for that loss.
It is not a calculation of expected persistence regret.

`n4_observability.py` enumerates multinomial count vectors for the full
record and binomial counts for a record reduced to K=4 versus K<4. It uses
floating-point finite sums, not simulation or a rational certificate.

| Years | Error using all hit counts | Error using only all-four years |
|---:|---:|---:|
| 25 | 0.3684 | 0.4781 |
| 50 | 0.3092 | 0.4725 |
| 100 | 0.2367 | 0.4597 |
| 200 | 0.1539 | 0.4423 |

The complete 100-year record is informative but leaves considerable error.
These numbers do not justify a claim that a century-long record cannot
distinguish the laws, or that m4 is practically unobservable. They also do
not show that real records suffice. Unknown lower moments, a continuous
class of candidate laws, observation error and dependent years change the
problem. A posterior for that problem requires an explicit prior and
observation model; none has been fitted here.

Checks cover the analytic one-year total-variation identity, identical-law
error of one half, relabelling count categories, normalization of both
enumerated likelihoods, and loss of information after reducing the record.
Reproduce with `uv run --python 3.12 --with numpy python ideation/b_screen/n4_observability.py`.

## Qualifications on the reversal map

The four fractions in `n4_reversal_map.json` count marked grid points:
389/1050, 9/1020, 264/1050 and 112/1020. They are not certified fractions
of continuous parameter-space area or probabilities of ecological
prevalence. Strict opposite signs imply local open reversal sets by
continuity, but an area bound needs their extents. Unmarked points do not
exclude reversals at other admissible laws. A causal explanation involving
background mortality remains a hypothesis, not a consequence of those
counts alone.

The exact n4 boundary-root brackets do support a third-moment ambiguity
band of width approximately 0.00127430303, about 8.191948056 percent of
the feasible m3 interval at the specified reversal parameters. This is a
fraction of interval length, not a probability mass over unknown laws.
