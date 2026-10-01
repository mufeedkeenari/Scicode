"""
Compute the factor multiplying the weight of the narrowest retained Gaussian of the 1/r^3 kernel splitting so that the short-range pressure kernel is continuous at the real-space cutoff.

The short-range part of the pressure kernel is truncated at a real-space cutoff r_cut. A
cutoff discontinuity in the instantaneous pressure produces artificial pressure jumps and
volume fluctuations in constant-pressure simulations, so the splitting is adjusted by
multiplying the weight of its narrowest retained Gaussian (index 0 of the retained terms)
by a dimensionless factor.

The continuity condition is imposed on the truncated long-range kernel built from exactly
the n_terms retained Gaussians, not on an untruncated series. The factor must be positive,
so that every long-range Gaussian keeps a positive weight.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_narrowest_weight_factor(b: float, sigma: float, n_terms: int, r_cut: float) -> float:
    '''Rescaling factor of the narrowest retained Gaussian for continuity at r_cut.

    Parameters
    ----------
    b : float
        Bilateral-series parameter b; must be finite and > 1.
    sigma : float
        Bilateral-series parameter sigma in Angstrom; must be finite and > 0.
    n_terms : int
        Number of retained long-range Gaussians; a positive integer.
    r_cut : float
        Real-space cutoff in Angstrom; must be finite and > 0.

    Returns
    -------
    omega : float
        Dimensionless factor multiplying the weight of the narrowest retained Gaussian.

    Raises
    ------
    ValueError
        If b, sigma or n_terms is invalid (as for the retained Gaussian terms), if r_cut
        is not a finite positive length, or if the continuity condition would require a
        non-positive factor.
    '''
    return omega

# =============================================================================
# GOLD SOLUTION
# =============================================================================

import numpy as np


def _oracle_compute_narrowest_weight_factor(b: float, sigma: float, n_terms: int, r_cut: float) -> float:
    if not (np.isfinite(r_cut) and r_cut > 0.0):
        raise ValueError("r_cut must be a finite positive length")
    weights, widths = _oracle_compute_pressure_kernel_gaussians(b, sigma, n_terms)
    tail = np.sum(weights[1:] * np.exp(-(r_cut / widths[1:]) ** 2))
    numerator = r_cut ** -3 - tail
    denominator = weights[0] * np.exp(-(r_cut / widths[0]) ** 2)
    if numerator <= 0.0 or denominator <= 0.0:
        raise ValueError("continuity at r_cut requires a non-positive narrowest weight")
    return float(numerator / denominator)

# =============================================================================
# TEST CASES
# =============================================================================

def test_cases():
    """Return list of test case specifications."""
    def invalid_case(args):
        setup = f"""def run_model():
    try:
        compute_narrowest_weight_factor({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_narrowest_weight_factor({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark splitting parameters.
        {
            "setup": "",
            "call": "compute_narrowest_weight_factor(1.6, 3.0, 10, 9.0)",
            "gold_call": "_oracle_compute_narrowest_weight_factor(1.6, 3.0, 10, 9.0)",
        },
        # Boundary: a single retained Gaussian carries the whole condition.
        {
            "setup": "",
            "call": "compute_narrowest_weight_factor(1.6, 3.0, 1, 9.0)",
            "gold_call": "_oracle_compute_narrowest_weight_factor(1.6, 3.0, 1, 9.0)",
        },
        # Edge: coarse ratio with only four retained terms, where the truncation of the
        # wide tail visibly shifts the factor.
        {
            "setup": "",
            "call": "compute_narrowest_weight_factor(2.0, 4.524309831775139, 4, 9.0)",
            "gold_call": "_oracle_compute_narrowest_weight_factor(2.0, 4.524309831775139, 4, 9.0)",
        },
        # Edge: shorter cutoff, factor well above one.
        {
            "setup": "",
            "call": "compute_narrowest_weight_factor(1.6, 3.2, 10, 8.0)",
            "gold_call": "_oracle_compute_narrowest_weight_factor(1.6, 3.2, 10, 8.0)",
        },
        # Invalid: non-positive cutoff.
        invalid_case("1.6, 3.0, 10, 0.0"),
        # Invalid: continuity would need a negative narrowest weight.
        invalid_case("2.0, 1.0, 30, 4.242640687119285"),
        # Invalid: upstream parameter check (b must exceed 1).
        invalid_case("0.5, 3.0, 10, 9.0"),
    ]
