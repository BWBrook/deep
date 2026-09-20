"""Regression test for persistence near numerical extinction."""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import b_screen_v2 as B


def test_long_horizon_endpoint_rankings_avoid_absorbing_mass_cancellation():
    """This cell was falsely classified as a reversal by ``1 - v[0]``."""
    p = 0.30
    q = p * p + (2 / 9) * (p - p * p)
    expected_deltas = (1.0811809108495736e-15, 1.0803951996431679e-15)

    for law, expected in zip((B.G_mix(0.0, p, q), B.G_mix(1.0, p, q)),
                             expected_deltas):
        _, persistence_r, persistence_c = B.run_case(
            4, 0.02, 0.95, 0.01, *law, T=200
        )
        delta = persistence_r - persistence_c
        assert delta > 0
        assert np.isclose(delta, expected, rtol=1e-12, atol=0.0)
