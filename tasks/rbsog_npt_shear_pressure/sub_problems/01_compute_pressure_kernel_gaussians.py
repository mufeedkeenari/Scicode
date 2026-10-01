"""
Compute the weights and widths of the retained long-range Gaussian terms in the bilateral-series sum-of-Gaussians approximation of the pressure-related 1/r^3 kernel.

The instantaneous Coulomb pressure tensor of a periodic charge system is governed by the
radial kernel 1/r^3. A smooth short-/long-range splitting of this kernel represents its
long-range part as a finite sum of Gaussians, each written as w_l * exp(-r^2 / s_l^2) with
weight w_l and width s_l (lengths in Angstrom, weights in Angstrom^-3). The splitting is
controlled by the bilateral-series parameters b > 1 and sigma > 0 (Angstrom) and by the
number n_terms of retained long-range Gaussians.

Return the retained terms ordered from the narrowest (index 0) to the widest
(index n_terms - 1), before any rescaling of the narrowest weight.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_pressure_kernel_gaussians(b: float, sigma: float, n_terms: int) -> "tuple[np.ndarray, np.ndarray]":
    '''Weights and widths of the retained long-range Gaussians of the 1/r^3 kernel.

    Parameters
    ----------
    b : float
        Bilateral-series parameter b; must be finite and > 1.
    sigma : float
        Bilateral-series parameter sigma in Angstrom; must be finite and > 0.
    n_terms : int
        Number of retained long-range Gaussians; a positive integer.

    Returns
    -------
    weights : np.ndarray
        Float array of shape (n_terms,), the weights w_l (before rescaling) in Angstrom^-3,
        ordered from the narrowest to the widest Gaussian.
    widths : np.ndarray
        Float array of shape (n_terms,), the widths s_l in Angstrom of the terms
        w_l * exp(-r^2 / s_l^2), in the same order.

    Raises
    ------
    ValueError
        If b is not a finite number greater than 1, sigma is not a finite positive
        number, or n_terms is not a positive integer.
    '''
    return weights, widths

# =============================================================================
# GOLD SOLUTION
# =============================================================================

import numpy as np


def _oracle_compute_pressure_kernel_gaussians(b: float, sigma: float, n_terms: int) -> "tuple[np.ndarray, np.ndarray]":
    if not (np.isfinite(b) and b > 1.0):
        raise ValueError("b must be a finite number greater than 1")
    if not (np.isfinite(sigma) and sigma > 0.0):
        raise ValueError("sigma must be a finite positive length")
    if isinstance(n_terms, bool) or not float(n_terms).is_integer() or n_terms < 1:
        raise ValueError("n_terms must be a positive integer")
    ell = np.arange(int(n_terms), dtype=float)
    weights = np.sqrt(2.0 / np.pi) * np.log(b) / (b ** (3.0 * ell) * sigma ** 3)
    widths = np.sqrt(2.0) * b ** ell * sigma
    return weights, widths

# =============================================================================
# TEST CASES
# =============================================================================

def test_cases():
    """Return list of test case specifications."""
    def invalid_case(args):
        setup = f"""def run_model():
    try:
        compute_pressure_kernel_gaussians({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_pressure_kernel_gaussians({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark splitting parameters (ten retained Gaussians).
        {
            "setup": "",
            "call": "compute_pressure_kernel_gaussians(1.6, 3.0, 10)",
            "gold_call": "_oracle_compute_pressure_kernel_gaussians(1.6, 3.0, 10)",
        },
        # Boundary: a single retained Gaussian.
        {
            "setup": "",
            "call": "compute_pressure_kernel_gaussians(1.6, 3.0, 1)",
            "gold_call": "_oracle_compute_pressure_kernel_gaussians(1.6, 3.0, 1)",
        },
        # Edge: coarse ratio b = 2 with a non-unit base length.
        {
            "setup": "",
            "call": "compute_pressure_kernel_gaussians(2.0, 4.524309831775139, 4)",
            "gold_call": "_oracle_compute_pressure_kernel_gaussians(2.0, 4.524309831775139, 4)",
        },
        # Edge: fine ratio close to 1, many terms.
        {
            "setup": "",
            "call": "compute_pressure_kernel_gaussians(1.1, 1.5, 40)",
            "gold_call": "_oracle_compute_pressure_kernel_gaussians(1.1, 1.5, 40)",
        },
        # Invalid: b must exceed 1.
        invalid_case("1.0, 3.0, 10"),
        # Invalid: sigma must be positive.
        invalid_case("1.6, 0.0, 10"),
        # Invalid: n_terms must be a positive integer.
        invalid_case("1.6, 3.0, 0"),
    ]
