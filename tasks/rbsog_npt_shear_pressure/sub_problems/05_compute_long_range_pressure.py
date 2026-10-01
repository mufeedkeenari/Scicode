"""
Compute the long-range (Fourier-space) contribution to the instantaneous Coulomb pressure tensor under the rescaled sum-of-Gaussians splitting of the 1/r^3 kernel, separated into its radial and non-radial parts.

The cell matrix has the lattice vectors as its columns and mode m labels the wavevector
k = 2 * pi * inv(cell).T @ m. The long-range part is the full sum over the nonzero integer
modes with max(|m_x|, |m_y|, |m_z|) <= m_max, with tinfoil boundary conditions (the m = 0
mode is omitted), for a charge-neutral system. Units follow the short-range step: charges in
e, lengths in Angstrom, unit Coulomb prefactor, pressure in e^2 / Angstrom^4. The long-range
Gaussians are the retained terms with the narrowest weight multiplied by its continuity
rescaling factor.

Each mode's contribution to the long-range tensor is a combination of the identity tensor
and the dyad k k. The radial part collects every contribution proportional to the identity
tensor and the non-radial part collects every contribution proportional to k k; their sum
is the full long-range pressure tensor.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_long_range_pressure(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int) -> "tuple[np.ndarray, np.ndarray]":
    '''Radial and non-radial parts of the long-range Coulomb pressure tensor.

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
        Mode truncation: all nonzero integer modes with max_d |m_d| <= m_max are summed;
        a positive integer.

    Returns
    -------
    radial : np.ndarray
        Float array of shape (3, 3), the radial part in e^2 / Angstrom^4 (a multiple of
        the identity tensor).
    nonradial : np.ndarray
        Symmetric float array of shape (3, 3), the non-radial part in e^2 / Angstrom^4.

    Raises
    ------
    ValueError
        If positions, charges or cell have invalid shapes or non-finite entries, if the
        cell determinant is not positive, if the system is not charge neutral, if m_max is
        not a positive integer, or if the splitting parameters are invalid (including a
        non-positive continuity rescaling factor).
    '''
    return radial, nonradial

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


def _validate_neutrality(charges):
    if abs(np.sum(charges)) > 1e-10 * max(1.0, np.sum(np.abs(charges))):
        raise ValueError("the system must be charge neutral")


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


def _oracle_compute_long_range_pressure(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int) -> "tuple[np.ndarray, np.ndarray]":
    positions, charges, cell = _validate_system(positions, charges, cell)
    _validate_neutrality(charges)
    modes, k, k2, spec5, spec7 = _mode_spectra(cell, b, sigma, n_terms, r_cut, m_max)
    power = _oracle_compute_structure_factor_power(positions, charges, cell, modes)
    volume = np.linalg.det(cell)
    radial = np.sum(spec5 * power) / (4.0 * volume ** 2) * np.eye(3)
    nonradial = -((spec7 * power)[:, None] * k).T @ k / (8.0 * volume ** 2)
    return radial, nonradial

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
        compute_long_range_pressure({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_long_range_pressure({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark snapshot with a converged mode set.
        {
            "setup": benchmark,
            "call": "compute_long_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6)",
            "gold_call": "_oracle_compute_long_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 6)",
        },
        # Boundary: smallest mode set (m_max = 1, 26 modes).
        {
            "setup": benchmark,
            "call": "compute_long_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 1)",
            "gold_call": "_oracle_compute_long_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0, 1)",
        },
        # Edge: dipolar ion pair in a strongly sheared cell with a narrow splitting, where
        # the non-radial tensor has large off-diagonal entries.
        {
            "setup": "import numpy as np\ncell = np.array([[8.0, 4.0, 2.0], [0.0, 7.0, 3.5], [0.0, 0.0, 9.0]])\npositions = np.array([[1.0, 1.0, 1.0], [4.0, 3.0, 5.5]])\ncharges = np.array([1.0, -1.0])\n",
            "call": "compute_long_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.5, 1.2, 12, 4.0, 8)",
            "gold_call": "_oracle_compute_long_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.5, 1.2, 12, 4.0, 8)",
        },
        # Invalid: net charge (tinfoil Fourier sum requires neutrality).
        invalid_case("np.zeros((2, 3)), np.array([1.0, 0.5]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 4"),
        # Invalid: non-positive mode truncation.
        invalid_case("np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]), np.array([1.0, -1.0]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 0"),
    ]
