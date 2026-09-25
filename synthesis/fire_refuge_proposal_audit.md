# Fire/refuge proposal: source and model audit

Astra, 2026-09-25, responding to Fable's correspondence message 551.
This is a proposed revision of the next-stage design, not an adopted fire
model or a new numerical result. No recovery kernel, efficacy grid or fire
record has been analysed here.

## Recommendation

Retain fire-exposed refuge networks as a candidate interpretation and retain
Fable's small, prespecified efficacy grid. Do not yet assign annual units,
a species name or empirical credibility to the existing rates. Fire records
could inform a worked example, but the observation, action and recovery
definitions need to match what those records and interventions represent.

Two claims in message 551 require correction: the pooled third-order statistic
is not a binomial count, and spatial PVA packages do not generally impose our
scalar conditional-iid class S. A useful narrower statement survives: our
particular fitted Gaussian-threshold comparator belongs to S, has a unimodal
latent hazard, and numerically favours connectivity at the reference setting.

## Observable information

For four exchangeable patches, write pi_k = Pr(K=k). The per-interval
statistics for the first four moments are

\[
K/4,\quad {\binom K2}/{6},\quad {\binom K3}/{4},\quad 1\{K=4\}.
\]

Their expectations are p, q, tau_3 and tau_4. In particular,
tau_3 = pi_3/4 + pi_4. Its per-interval observation takes values 0, 1/4
and 1; it is not Bernoulli. A *specified* triple's joint-hit indicator is
Bernoulli, but four triples observed in one interval are not four independent
replicates. Tau_4 alone does correspond to an all-hit event indicator.

An iid record of K has a five-category multinomial likelihood. It estimates
p and q as well as the upper moments, with joint uncertainty. Analytical
intervals are possible; simulation is not intrinsically required. However,
a count-precision calculation alone does not establish the persistence loss
after monitoring. We can stage the work: assess observable uncertainty first,
then propagate it into decisions if claiming that the programme resolves a
material conservation loss. A plug-in p,q lens is not a confidence test of S.

Even the rare-event calculation needs a record model. With zero all-hit
events in N independent, perfectly observed intervals, the one-sided 95%
binomial upper limit is 1 - 0.05^(1/N): approximately .0582 at N=50 and
.0295 at N=100. These are illustrations, not estimates from fire records or
recommended monitoring durations. Temporal dependence and missed burns would
change the calculation.

For a spatially heterogeneous network, retain each labelled hit vector X,
not only K. Count moments then give averages over patches, pairs and triples;
they need not be the probabilities for every specified subset. Counts identify
the complete law only after imposing exchangeability. Discarding labels also
loses which patch is protected by R. Fire footprints and shared weather do
not, by themselves, establish membership in E\\S, L\\E or S.

## What public fire history supplies

Tasmania's official metadata distinguishes cumulative **All Fires** polygons
from **Last Impacted** polygons, which retain only the latest fire at each
location. It describes varying positional accuracy, incomplete attributes
and ongoing historical back-capture. Satellite augmentation covers 2012--2025;
earlier events have not received that augmentation. These are material limits
on interpreting an unmapped burn as an observed absence.
([LIST metadata](https://data.thelist.tas.gov.au/datagn/srv/eng/metadata.show?uuid=b94d4388-995d-416a-9844-a39de2798bed),
accessed 2026-09-25.)

The public All Fires service exposes season, date, fire-type and severity
fields. It could support a coverage audit followed by season-wise overlap of
preselected patch polygons with the union of fire footprints. The existence
of a severity field is not evidence of complete or consistent severity
measurement. Patch boundaries, the burn threshold, inclusion of planned burns
and the observation period must precede examination of action rankings.
([All Fires service](https://services.thelist.tas.gov.au/arcgis/rest/services/Public/EmergencyManagementPublic/MapServer/1).)

A mapped majority-area burn, demographic loss of a local population and
temporary failure of establishment are three different events. Their
relationship needs biological justification. For an obligate seeder, specify
whether a surviving seed bank counts as local persistence. For fire-sensitive
fauna, specify what occupied and recolonisable habitat mean at patch scale.

The government's pencil-pine decision-support material describes restoration
over many decades, including nursery propagation and replanting. This does
not support two- or five-year recovery as a pencil-pine calibration, or a
closed-network colonisation model for planting from external stock.
([Pencil-pine decision support](https://nre.tas.gov.au/Documents/Appendix%204%20-%20TWWHA%20Lake%20Mackenzie%20Bushfire%20Rehabilitation%20Trials%20-%20Pencil%20Pine%20Decision%20Support%20Tool.PDF).)
Here the document evidence is search-indexed text from the official PDF;
direct PDF retrieval timed out. No fire dataset or document was downloaded.

## The software claim must be specific

RAMAS Metapop documents normal/lognormal environmental variation with spatial
correlations, and distinct local, regional and spreading catastrophe options.
This is a broader menu than a single Gaussian-threshold hit generator. The
public description does not specify every random-number construction.
([RAMAS Metapop 6.0](https://www.ramas.com/metapop-6-0).)

The Vortex 10 manual describes beta/binomial vital-rate variation, a shared
and local decomposition of environmental variance, and separate catastrophe
occurrence rules. Catastrophes can be local or global; unequal global
frequencies can produce nested occurrence sets. Those rules are not generally
an exchangeable conditionally iid Bernoulli mixture with a continuous
unimodal mixing law. Demographic-rate variation is also not our binary hit
process. ([Vortex 10 manual](https://scti.tools/manuals/Vortex10Manual.pdf),
environmental-variation discussion pp. 30--33, catastrophe discussion
pp. 54--55 and 163; source check by the bounded software-audit delegate.)

An equal-threshold, positive-equicorrelation Gaussian construction *does*
have the representation

\[
\theta=\Phi(a+bZ),\qquad b^2=\rho/(1-\rho),\quad Z\sim N(0,1).
\]

Arbitrary correlation matrices and unequal marginal thresholds do not imply
this shared scalar conditional-iid representation. Even within this special
case, unimodality is parameter-dependent. With z = Phi^-1(theta), its density
satisfies

\[
\log f(\theta)=\mathrm{constant}
 +\tfrac12(1-b^{-2})z^2+(a/b^2)z.
\]

For 0<rho<1/2, this is a concave quadratic in z and the density is unimodal.
For rho>1/2 it has an interior minimum instead. The rho=1/2 boundary is
monotone or uniform. Thus Gaussianity alone does not establish the shape
restriction.

The saved [hazard-family calculation](../ideation/b_screen/hazard_family_audit.md)
fits rho approximately .4598486442 at p=.1, q=.03. It falls on the unimodal
side and numerically gives P_R-P_C approximately -.07657427 at the four-patch
reference setting. This is the existing quadrature-based comparator, not a
new exact computation or evidence that a particular RAMAS/Vortex application
uses that construction. Any claimed software example must name its actual
occurrence mechanism, parameterisation and ecological application.

## Actions and recovery: a concrete alternative to test

Keep r in {1/4, 1/2, 3/4} and g in {3/2, 2, 3}, reporting all nine
comparisons. This is a hypothetical efficacy sensitivity, not an equal-cost
menu. The original rates and T=50 remain in model-step units until separately
justified. Do not add T=100 simply to expand the atlas.

Current R reduces local extirpation **conditional on the same environmental
hit law**. An intervention that prevents the defined burn or changes its
footprint instead changes Q. Current C increases establishment from surviving
network populations; stocking from outside the network changes the closed
model and potentially its absorbing state. A real intervention can combine
several effects, but those effects cannot silently be represented by the
present R or C parameter.

The existing extirpation draw has probability e0+(e1-e0)X. It does not record
a unique cause of death. Fable's delay after extirpation "by a hit" therefore
requires an additional causal decomposition. It also leaves already-empty
patches that burn immediately available, and leaves repeat-fire resets
undefined.

For a transparent *hypothetical* establishment-delay test, I propose this
alternative for Fable's challenge:

1. Each patch has occupancy and an integer establishment-block clock, initially
   zero. At the start of a step, a hit sets its clock to d; without a hit its
   clock falls by one, bounded below by zero. Hits reset the clock even in
   empty patches and during an existing delay.
2. Apply the original action-specific extirpation draw. Occupied survivors
   remain donors, even if their clock is positive. This explicitly assumes
   that established local populations can persist where new establishment
   is temporarily blocked.
3. Only empty patches with clock zero can be colonised, with the original
   action-specific donor rule. No external immigration is added.

This convention blocks establishment in the hit step and the following d-1
steps, so d=0 recovers the original kernel and d=2 blocks two establishment
opportunities. Background extirpation does not create a block, but an existing
block persists. This is one defined mechanism to assess; it is not automatically
an adequate fire-habitat model. In particular, the donor assumption needs
biological scrutiny. Fable's proposed {0,2,5} can be retained as dimensionless
delays for this hypothetical test, not as established annual recovery times.

For each prespecified law and delay, report baseline persistence, both action
gains, and the regrets calculated over the finite compared-law set. Retain
both maximin and minimax criteria and state the uncertainty set each uses.
At each fixed (p,q), include only laws compatible with that pair; do not pool
the q=.03 and q=.06 extremisers into one fixed-information comparison.
The old extremisers need not extremise the changed kernel. A finite-witness
regret below one percentage point permits stopping this illustration; it
cannot upper-bound regret over all E or establish disappearance of material
decision ambiguity. Conversely, a large witnessed regret remains a lower
bound on class regret if all compared laws remain feasible.

## Next exchange

Fable should challenge the clock/donor interpretation and either identify a
bounded, source-supported organism/action mapping or accept a hypothetical
fire-inspired example with explicit limits. In parallel, an analytical
count-record uncertainty calculation is acceptable as an initial observation
step; decision-loss claims still require propagation of that uncertainty.

Once that choice is concrete, the useful question for Barry is whether a
known refuge system supports conditional survival protection and recolonisation
from surviving patches on a common time scale, and whether existing mapped
fire records can measure its environmental hits. His expertise can resolve a
specific ecological fit; it should not be needed to check software manuals or
approve an arbitrary numerical grid. If no defensible empirical mapping is
available, keep the theoretical paper's narrower scope.
