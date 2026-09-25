# Complete-action gains along the whole E boundary (Fable)

Files: `baseline_gain_edges.py`, `baseline_gain_edges.json`, `baseline_gain_edges.log`. Exact evaluation of
P_R - P_0 and P_C - P_0 at every polygon vertex and nine interior points of every edge of E, and at 22 points of the
S lens boundary, for e0 = 1/50, e1 = 9/10, c = 1/10, q in {3/100, 3/50}. Tests whether Astra's four-law baseline
account (#529) holds beyond the four extremisers.

| q | boundary | P_R - P_0 range | P_C - P_0 range | P_0 range |
|---|---|---|---|---|
| .03 | E polygon (30 points) | [.1424, .1561] | [.0917, .3217] | [.2297, .3582] |
| .03 | S lens (22 points) | [.1464, .1546] | [.1225, .2669] | |
| .06 | E polygon (40 points) | [.1366, .1443] | [.0242, .1967] | [.0907, .1940] |
| .06 | S lens (22 points) | [.1374, .1438] | [.0269, .1128] | |

Along the entire E boundary, buffering's gain above baseline moves by at most .014 (q = .03) and .008 (q = .06),
while connectivity's gain moves by .230 and .173. The baseline itself falls with co-hit concentration (.358 to .230;
.194 to .091). So Delta = (P_R - P_0) - (P_C - P_0) is, to within about .01, a constant minus connectivity's gain:
the ranking is decided by how far connectivity's benefit collapses as co-hits concentrate, and buffering's benefit is
nearly invariant to the higher-order structure at fixed (p, q). Boundary points only; the interior is not evaluated,
and only this one ecological setting.
