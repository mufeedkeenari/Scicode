"""
Compute the single-mode variance of the random-batch estimate of one Cartesian component of the non-radial long-range Coulomb pressure tensor when the Fourier mode is drawn from the non-radial importance distribution.

A random-batch estimate of the non-radial long-range pressure averages an unbiased
single-mode estimate over P independently drawn modes, so its variance is the single-mode
variance divided by P. This step returns that single-mode (P = 1) variance for the (mu, nu)
component, with the mode drawn exactly from the non-radial importance distribution on the
nonzero integer modes with max(|m_x|, |m_y|, |m_z|) <= m_max (mode m labels the wavevector
k = 2 * pi * inv(cell).T @ m). The estimate's mean is the corresponding component of the
non-radial long-range tensor summed over the same mode set.

Units: charges in e, lengths in Angstrom, unit Coulomb prefactor; the variance is in
e^4 / Angstrom^8. Cartesian indices are (x, y, z) = (0, 1, 2).
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_nonradial_variance(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int, mu: int, nu: int) -> float:
    '''Single-mode variance of the non-radial random-batch pressure estimate.

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
        Real-space cutoff in Angstrom used for the continuity rescaling; finite and > 0.
    m_max : int
        Mode truncation (nonzero modes with max_d |m_d| <= m_max); a positive integer.
    mu : int
        First Cartesian index of the tensor component, in {0, 1, 2}.
    nu : int
        Second Cartesian index of the tensor component, in {0, 1, 2}.

    Returns
    -------
    variance : float
        Single-mode variance of the (mu, nu) estimate in e^4 / Angstrom^8.

    Raises
    ------
    ValueError
        If mu or nu is not an integer in {0, 1, 2}, if positions, charges or cell are
        invalid, if the system is not charge neutral, if m_max is not a positive integer,
        or if the splitting parameters are invalid (including a non-positive continuity
        rescaling factor).
    '''
    return variance

# =============================================================================
# GOLD SOLUTION
# =============================================================================

import numpy as np


def _validate_cell(cell):
    cell = np.array(cell, dtype=float)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
        raise ValueError("cell must be a finite 3x3 matrix")
    if np.linalg.det(cell) <= 0.0:
        raise ValueError("cell must have a positive determinant")
    return cell


def _validate_system(positions, charges, cell):
    positions = np.array(positions, dtype=float)
    charges = np.array(charges, dtype=float)
    cell = _validate_cell(cell)
    if positions.ndim != 2 or positions.shape[1] != 3 or positions.shape[0] < 1:
        raise ValueError("positions must have shape (N, 3) with N >= 1")
    if charges.shape != (positions.shape[0],):
        raise ValueError("charges must have shape (N,)")
    if not (np.all(np.isfinite(positions)) and np.all(np.isfinite(charges))):
        raise ValueError("positions and charges must be finite")
    return positions, charges, cell


def _rescaled_gaussians(b, sigma, n_terms, r_cut):
    weights, widths = _oracle_compute_pressure_kernel_gaussians(b, sigma, n_terms)
    weights = weights.copy()
    weights[0] *= _oracle_compute_narrowest_weight_factor(b, sigma, n_terms, r_cut)
    return weights, widths


def _integer_modes(m_max):
    if isinstance(m_max, bool) or not float(m_max).is_integer() or m_max < 1:
        raise ValueError("m_max must be a positive integer")
    axis = np.arange(-int(m_max), int(m_max) + 1)
    modes = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    return modes[np.any(modes != 0, axis=1)]


def _mode_spectra(cell, b, sigma, n_terms, r_cut, m_max):
    """Modes, wavevectors, |k|^2 and the s^5 / s^7 Gaussian spectra of the kernel."""
    modes = _integer_modes(m_max)
    weights, widths = _rescaled_gaussians(b, sigma, n_terms, r_cut)
    k = 2.0 * np.pi * modes @ np.linalg.inv(cell)
    k2 = np.sum(k ** 2, axis=1)
    gauss = np.exp(-0.25 * k2[:, None] * widths ** 2)
    spec5 = np.pi ** 1.5 * gauss @ (weights * widths ** 5)
    spec7 = np.pi ** 1.5 * gauss @ (weights * widths ** 7)
    return modes, k, k2, spec5, spec7


def _oracle_compute_nonradial_variance(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int, mu: int, nu: int) -> float:
    for index in (mu, nu):
        if isinstance(index, bool) or not float(index).is_integer() or index not in (0, 1, 2):
            raise ValueError("mu and nu must be Cartesian indices in {0, 1, 2}")
    mu, nu = int(mu), int(nu)
    positions, charges, cell = _validate_system(positions, charges, cell)
    _, nonradial = _oracle_compute_long_range_pressure(positions, charges, cell, b, sigma, n_terms, r_cut, m_max)
    normalization = _oracle_compute_nonradial_normalization(cell, b, sigma, n_terms, r_cut, m_max)
    modes, k, k2, _, spec7 = _mode_spectra(cell, b, sigma, n_terms, r_cut, m_max)
    power = _oracle_compute_structure_factor_power(positions, charges, cell, modes)
    volume = np.linalg.det(cell)
    probability = spec7 * k2 ** 2 / normalization
    sample = -normalization * power * k[:, mu] * k[:, nu] / (8.0 * volume ** 2 * k2 ** 2)
    second_moment = np.sum(probability * sample ** 2)
    return float(max(second_moment - nonradial[mu, nu] ** 2, 0.0))

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
        compute_nonradial_variance({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_nonradial_variance({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark shear component (x, y).
        {
            "setup": benchmark,
            "call": "compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 0, 1)",
            "gold_call": "_oracle_compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 0, 1)",
        },
        # Edge: diagonal component (z, z) of the same snapshot.
        {
            "setup": benchmark,
            "call": "compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 2, 2)",
            "gold_call": "_oracle_compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6, 2, 2)",
        },
        # Boundary: smallest mode set (m_max = 1) for the (y, z) component.
        {
            "setup": benchmark,
            "call": "compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 1, 1, 2)",
            "gold_call": "_oracle_compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 1, 1, 2)",
        },
        # Edge: dipolar ion pair in a strongly sheared cell, (x, z) component.
        {
            "setup": "import numpy as np\ncell = np.array([[8.0, 4.0, 2.0], [0.0, 7.0, 3.5], [0.0, 0.0, 9.0]])\npositions = np.array([[1.0, 1.0, 1.0], [4.0, 3.0, 5.5]])\ncharges = np.array([1.0, -1.0])\n",
            "call": "compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.5, 1.2, 12, 4.0, 8, 0, 2)",
            "gold_call": "_oracle_compute_nonradial_variance(positions.copy(), charges.copy(), cell.copy(), 1.5, 1.2, 12, 4.0, 8, 0, 2)",
        },
        # Invalid: Cartesian index out of range.
        invalid_case("np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]), np.array([1.0, -1.0]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 4, 0, 3"),
        # Invalid: net charge.
        invalid_case("np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]), np.array([1.0, 1.0]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 4, 0, 1"),
    ]
