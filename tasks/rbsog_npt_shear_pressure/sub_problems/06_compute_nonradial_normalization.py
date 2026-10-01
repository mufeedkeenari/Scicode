"""
Compute the normalization constant of the non-radial importance-sampling distribution over Fourier modes used by the random-batch estimator of the non-radial long-range pressure.

The random-batch method replaces the full Fourier sum of the non-radial long-range pressure
by an average over a small batch of modes drawn from a discrete importance distribution
defined by the rescaled long-range Gaussians. The distribution is supported on the
nonzero integer modes m with max(|m_x|, |m_y|, |m_z|) <= m_max, where mode m labels the
wavevector k = 2 * pi * inv(cell).T @ m for a cell matrix whose columns are the lattice
vectors. The returned constant is the normalization of that distribution, summed over
exactly this truncated mode set. The importance distribution depends only on the cell and
the splitting, not on the particle configuration.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_nonradial_normalization(cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int) -> float:
    '''Normalization constant of the non-radial Fourier-mode importance distribution.

    Parameters
    ----------
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
        Mode truncation: the distribution is supported on all nonzero integer modes with
        max_d |m_d| <= m_max; a positive integer.

    Returns
    -------
    normalization : float
        The (dimensionless) normalization constant of the non-radial distribution.

    Raises
    ------
    ValueError
        If cell is not a finite 3x3 matrix with positive determinant, if m_max is not a
        positive integer, or if the splitting parameters are invalid (including a
        non-positive continuity rescaling factor).
    '''
    return normalization

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


def _oracle_compute_nonradial_normalization(cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float, m_max: int) -> float:
    cell = _validate_cell(cell)
    _, _, k2, _, spec7 = _mode_spectra(cell, b, sigma, n_terms, r_cut, m_max)
    return float(np.sum(spec7 * k2 ** 2))

# =============================================================================
# TEST CASES
# =============================================================================

def test_cases():
    """Return list of test case specifications."""
    def invalid_case(args):
        setup = f"""import numpy as np
def run_model():
    try:
        compute_nonradial_normalization({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_nonradial_normalization({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    cell = "import numpy as np\ncell = np.array([[9.6, 1.4, -0.7], [0.0, 9.1, 1.1], [0.0, 0.0, 10.2]])\n"
    return [
        # Normal: benchmark cell and splitting with a converged mode set.
        {
            "setup": cell,
            "call": "compute_nonradial_normalization(cell.copy(), 1.6, 3.0, 10, 9.0, 6)",
            "gold_call": "_oracle_compute_nonradial_normalization(cell.copy(), 1.6, 3.0, 10, 9.0, 6)",
        },
        # Boundary: smallest mode set (m_max = 1).
        {
            "setup": cell,
            "call": "compute_nonradial_normalization(cell.copy(), 1.6, 3.0, 10, 9.0, 1)",
            "gold_call": "_oracle_compute_nonradial_normalization(cell.copy(), 1.6, 3.0, 10, 9.0, 1)",
        },
        # Edge: large cubic cell with a narrow splitting, where many modes contribute and
        # the wide Gaussians are resolved by the reciprocal lattice.
        {
            "setup": "import numpy as np\ncell = 30.0 * np.eye(3)\n",
            "call": "compute_nonradial_normalization(cell.copy(), 1.3, 1.0, 20, 3.5, 12)",
            "gold_call": "_oracle_compute_nonradial_normalization(cell.copy(), 1.3, 1.0, 20, 3.5, 12)",
        },
        # Invalid: singular cell.
        invalid_case("np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 0.0, 1.0]]), 1.6, 3.0, 10, 9.0, 4"),
        # Invalid: non-integer mode truncation.
        invalid_case("10.0 * np.eye(3), 1.6, 3.0, 10, 9.0, 2.5"),
    ]
