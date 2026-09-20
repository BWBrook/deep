# Full-text antecedent audit, Astra, 2026-09-20

This note records an independent reading of the publisher PDFs supplied by
Barry, filed by Fable, and listed in `resources/LEDGER.md`. Both PDF hashes
match the ledger. Text extraction was checked against rendered pages for
the principal passages below. Derived text and page images are untracked
under `resources/literature/`; original PDFs were not changed.

## Ferson and Burgman (1995), Biological Conservation 73:101-105

DOI: https://doi.org/10.1016/0006-3207(95)90031-4

The paper is a direct conceptual ancestor. Page 103 explicitly distinguishes
joint distributions that share means, standard deviations, and correlations;
it then bounds risk under unspecified dependence. Pages 103-104 show why
extreme linear-correlation assumptions need not enclose the risk arising from
more general dependence. Calling this "no overlap with the identification
framing" (Fable, correspondence #372) is too strong.

The actual worked dependency bounds fix marginal vital-rate distributions
and leave their joint law unrestricted. The example is a stage-structured
single population, with one-year quasi-extinction risk as its main target.
The five-year temporal-dependence comparison is partial: the paper explicitly
explains the difficulty of obtaining optimal dependency bounds after repeated
appearances of the same random quantities. It does not provide our specified
fixed-bivariate-binary-law moment problem, finite-horizon patch Markov model,
or buffering-versus-connectivity ranking certificate.

The candidate contribution should therefore be described as a specific
extension and sharpening of dependency-bound reasoning, not the invention
of an identification viewpoint on ecological risk. These distinctions are
supported by this paper; they do not establish absence of an exact precedent
elsewhere in the literature.

## Heino, Kaitala, Ranta and Lindstrom (1997), Proc. R. Soc. B 264:481-486

DOI: https://doi.org/10.1098/rspb.1997.0069

The paper combines Ricker abundance dynamics, distance-dependent deterministic
migration, local noise, and an intermittent common disturbance. Its synchrony
measure is average pairwise correlation of abundance time series (methods,
pp. 481-482), not a marginal/pairwise binary environmental-hit law.
Populations are reinitialized after global extinction, so its repeated-event
extinction rate is also distinct from our absorbing-chain finite-horizon
survival probability.

The common-driver component is related to our motivation, but the complete
population process is not our conditional-iid homogeneous mortality model:
patch states, nonlinear growth, and spatial dispersal affect subsequent
outcomes. Conditioning only on the common disturbance does not remove these
features or make the local population dynamics identical.

The discussion (pp. 485-486) explains why its deterministic migration enables
recolonisation whenever a population survives and why migration can both
rescue and synchronise populations. Our model holds the environmental hazard
law fixed when increasing connectivity, and contains no feedback from
connectivity to environmental synchrony. Thus its management implications
cannot simply be transferred between these models.

The reported close association between synchrony and extinction under the
specified simulations is useful context. It is not a theorem that pairwise
synchrony identifies regional risk whenever a mechanism is known or described
by one parameter; that requires an invertible mapping from the measured
statistic to the remaining relevant model quantities.

## Status

These two full-text gaps are closed. The work remains motivated by established
dependency and synchrony research. Novelty of the specific fixed-pairwise-law
action problem remains a candidate claim, to be tested through the broader
antecedent search rather than inferred from model agreement.
