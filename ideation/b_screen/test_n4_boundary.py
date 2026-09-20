"""Focused regression tests for exact boundary certificates."""
from fractions import Fraction as F
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import n4_boundary as B


def test_depth_zero_preserves_valid_outer_bounds_and_reports_unmet_tolerance():
    # Bernstein coefficients [0, -1, 0] represent -2*u*(1-u), whose interior
    # minimum is -1/2 while endpoint evaluation alone attains 0.
    result = B.enclose([B.Q(0), B.Q(-1), B.Q(0)], tol=F(1, 10**9), maxd=0)
    minimum, maximum = result["min"], result["max"]
    assert minimum["outer"] == B.Q(-1)
    assert minimum["attained"] == B.Q(0)
    assert minimum["witness_u"] == F(0)
    assert not minimum["tolerance_met"]
    assert minimum["depth_limit_reached"]
    assert maximum["outer"] == B.Q(0)
    assert maximum["attained"] == B.Q(0)


def test_converged_certificate_has_exact_interior_witness_and_valid_interval():
    result = B.enclose([B.Q(0), B.Q(-1), B.Q(0)], tol=F(1, 10**9), maxd=20)
    minimum, maximum = result["min"], result["max"]
    assert minimum["attained"] == B.Q(-1, 2)
    assert minimum["witness_u"] == F(1, 2)
    assert minimum["outer"] <= B.Q(-1, 2) <= minimum["attained"]
    assert maximum["outer"] >= B.Q(0) >= maximum["attained"]
    assert minimum["tolerance_met"] and maximum["tolerance_met"]


def test_float_compatibility_endpoints_are_outward():
    exact = B.Q(1, 10)
    down = B._float_outward(exact, "down")
    up = B._float_outward(exact, "up")
    assert B.Q(down) <= exact
    assert B.Q(up) >= exact
    assert down == math.nextafter(float(exact), float("-inf"))
    assert up == float(exact)


def test_straddling_bounds_with_same_sign_attained_values_are_unresolved():
    minimum = {"outer": B.Q(-1, 10), "attained": B.Q(1, 10)}
    maximum = {"outer": B.Q(3, 10), "attained": B.Q(2, 10)}
    assert B.classify_whole_class(minimum, maximum) == "UNRESOLVED"
    minimum["attained"] = B.Q(-1, 10)
    assert B.classify_whole_class(minimum, maximum) == "UNIDENTIFIED"

