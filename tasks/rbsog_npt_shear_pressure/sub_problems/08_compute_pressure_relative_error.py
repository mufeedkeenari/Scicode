"""
Compute the relative standard error of the random-batch sum-of-Gaussians estimate of an off-diagonal (shear) component of the instantaneous Coulomb pressure tensor for a given batch size (end-to-end pipeline; orchestrator).

The deterministic reference value of the (mu, nu) shear component is the rescaled
sum-of-Gaussians pressure tensor: its short-range real-space part plus its full long-range
Fourier-space part (radial plus non-radial) on the truncated mode set. In the random-batch
method the short-range part is evaluated exactly, while the long-range part is replaced by
importance-sampled mini-batches of Fourier modes. Here the batch_size modes of the
non-radial estimate are treated as independent draws from the non-radial importance
distribution, i.e. an ideally mixed measure-recalibration chain.

The returned quantity is the standard deviation of the batch estimate of the (mu, nu)
component divided by the magnitude of its deterministic reference value (dimensionless).
Only off-diagonal components (mu != nu) are accepted.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_pressure_relative_error(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int, batch_size: int, mu: int, nu: int) -> float:
    '''Relative standard error of the random-batch estimate of a shear pressure component.

    Parameters
    ----------
    positions : np.ndarray
        Cartesian positions in Angstrom, shape (N, 3) with N >= 1. Not modified.
    charges : np.ndarray
        Charges in units of e, shape (N,); the total charge must vanish (to within
        1e-10 times max(1, sum |q_i|)). Not modified.
    cell : np.ndarray
        Cell matrix of shape (3, 3) whose columns are the lattice vectors in Angstrom;
        its determinant must be positive. Not modified.
    b : float
        Bilateral-series parameter b; must be finite and > 1.
    sigma : float
        Bilateral-series parameter sigma in Angstrom; must be finite and > 0.
    n_terms : int
        Number of retained long-range Gaussians; a positive integer.
    r_cut : float
        Real-space cutoff in Angstrom; must be finite and > 0.
    m_max : int
        Mode truncation (nonzero modes with max_d |m_d| <= m_max); a positive integer.
    batch_size : int
        Number P of independently drawn non-radial modes; a positive integer.
    mu : int
        First Cartesian index of the shear component, in {0, 1, 2}.
    nu : int
        Second Cartesian index of the shear component, in {0, 1, 2}, with nu != mu.

    Returns
    -------
    relative_error : float
        Standard deviation of the batch estimate of the (mu, nu) pressure component
        divided by the magnitude of its deterministic sum-of-Gaussians value.

    Raises
    ------
    ValueError
        If mu == nu, if batch_size is not a positive integer, if the deterministic
        (mu, nu) component is exactly zero, or for any invalid input rejected by the
        short-range, long-range or variance computations (including a non-neutral
        system and invalid indices).
    '''
    return relative_error

# =============================================================================
# GOLD SOLUTION
# =============================================================================

import numpy as np


def _oracle_compute_pressure_relative_error(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int, batch_size: int, mu: int, nu: int) -> float:
    if mu == nu:
        raise ValueError("an off-diagonal component (mu != nu) is required")
    if isinstance(batch_size, bool) or not float(batch_size).is_integer() or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    variance = _oracle_compute_nonradial_variance(positions, charges, cell, b, sigma, n_terms, r_cut, m_max, mu, nu)
    short = _oracle_compute_short_range_pressure(positions, charges, cell, b, sigma, n_terms, r_cut)
    radial, nonradial = _oracle_compute_long_range_pressure(positions, charges, cell, b, sigma, n_terms, r_cut, m_max)
    reference = short[mu, nu] + radial[mu, nu] + nonradial[mu, nu]
    if reference == 0.0:
        raise ValueError("the deterministic pressure component vanishes")
    return float(np.sqrt(variance / batch_size) / abs(reference))

# =============================================================================
# TEST CASES
# =============================================================================

def test_cases():
    """Return list of test case specifications."""
    benchmark = """import numpy as np
cell = np.array([[9.6, 1.4, -0.7], [0.0, 9.1, 1.1], [0.0, 0.0, 10.2]])
frac = np.array([[0.02, 0.03, 0.01], [0.47, 0.55, 0.04], [0.53, 0.06, 0.46], [0.05, 0.49, 0.52],
                 [0.51, 0.02, 0.07], [0.03, 0.54, 0.03], [0.06, 0.03, 0.55], [0.48, 0.46, 0.51]])
positions = frac @ cell.T
charges = np.array([1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0])
"""
    def invalid_case(args):
        setup = f"""import numpy as np
def run_model():
    try:
        compute_pressure_relative_error({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_pressure_relative_error({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark shear component P_xy at batch size 128.
        {
            "setup": benchmark,
            "call": "compute_pressure_relative_error(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)",
            "gold_call": "_oracle_compute_pressure_relative_error(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)",
        },
        # Boundary: a single sampled mode (batch_size = 1), component (z, x).
        {
            "setup": benchmark,
            "call": "compute_pressure_relative_error(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 1, 2, 0)",
            "gold_call": "_oracle_compute_pressure_relative_error(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 1, 2, 0)",
        },
        # Edge: dipolar ion pair in a strongly sheared cell with a short cutoff.
        {
            "setup": "import numpy as np\ncell = np.array([[8.0, 4.0, 2.0], [0.0, 7.0, 3.5], [0.0, 0.0, 9.0]])\npositions = np.array([[1.0, 1.0, 1.0], [4.0, 3.0, 5.5]])\ncharges = np.array([1.0, -1.0])\n",
            "call": "compute_pressure_relative_error(positions.copy(), charges.copy(), cell.copy(), 1.5, 1.2, 12, 4.0, 8, 64, 0, 2)",
            "gold_call": "_oracle_compute_pressure_relative_error(positions.copy(), charges.copy(), cell.copy(), 1.5, 1.2, 12, 4.0, 8, 64, 0, 2)",
        },
        # Invalid: diagonal component.
        invalid_case("np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]), np.array([1.0, -1.0]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 4, 16, 1, 1"),
        # Invalid: non-positive batch size.
        invalid_case("np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]), np.array([1.0, -1.0]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 4, 0, 0, 1"),
    ]
