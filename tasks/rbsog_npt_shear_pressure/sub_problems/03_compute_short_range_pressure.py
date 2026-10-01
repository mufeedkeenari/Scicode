"""
Compute the short-range (real-space) contribution to the instantaneous Coulomb pressure tensor of a periodic point-charge system under the rescaled sum-of-Gaussians splitting of the 1/r^3 kernel.

The system is a triclinic periodic cell whose 3x3 cell matrix has the lattice vectors as
its columns, so a particle with fractional coordinates s sits at r = cell @ s and lattice
translations are cell @ n for integer vectors n. Charges are in units of e and lengths in
Angstrom, with Gaussian electrostatic units and a unit Coulomb prefactor, so the pressure
tensor is returned in e^2 / Angstrom^4.

The short-range kernel is truncated at the real-space cutoff r_cut and uses the retained
long-range Gaussians with the narrowest weight multiplied by its continuity rescaling
factor. Interactions are summed over every periodic image separated by less than r_cut,
including images of a particle with itself; the cutoff may exceed half the cell width, so
a minimum-image sum is not sufficient. Only a particle's interaction with itself in the
same image is excluded.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_short_range_pressure(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float) -> "np.ndarray":
    '''Short-range Coulomb pressure tensor of the rescaled SOG splitting.

    Parameters
    ----------
    positions : np.ndarray
        Cartesian positions in Angstrom, shape (N, 3) with N >= 1; positions need not lie
        inside the reference cell. Not modified.
    charges : np.ndarray
        Charges in units of e, shape (N,). Not modified.
    cell : np.ndarray
        Cell matrix of shape (3, 3) whose columns are the lattice vectors in Angstrom;
        its determinant (the cell volume) must be positive. Not modified.
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
    pressure : np.ndarray
        Symmetric float array of shape (3, 3), the short-range pressure tensor in
        e^2 / Angstrom^4, indexed by Cartesian components (x, y, z) = (0, 1, 2).

    Raises
    ------
    ValueError
        If positions, charges or cell have invalid shapes or non-finite entries, if the
        cell determinant is not positive, or if the splitting parameters are invalid
        (including a non-positive continuity rescaling factor).
    '''
    return pressure

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


def _oracle_compute_short_range_pressure(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", b: float, sigma: float, n_terms: int, r_cut: float) -> "np.ndarray":
    positions, charges, cell = _validate_system(positions, charges, cell)
    weights, widths = _rescaled_gaussians(b, sigma, n_terms, r_cut)
    inv_cell = np.linalg.inv(cell)
    frac = positions @ inv_cell.T
    dfrac = frac[:, None, :] - frac[None, :, :]
    dfrac -= np.round(dfrac)
    # Image range covering every separation shorter than r_cut.
    reach = np.ceil(r_cut * np.linalg.norm(inv_cell, axis=1)).astype(int) + 1
    axes = [np.arange(-n, n + 1) for n in reach]
    images = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    same_image = np.all(images == 0, axis=1)
    pressure = np.zeros((3, 3))
    for i in range(len(charges)):
        sep = (dfrac[i][:, None, :] + images[None, :, :]) @ cell.T
        dist = np.linalg.norm(sep, axis=-1)
        inside = dist < r_cut
        inside[i] &= ~same_image
        d = sep[inside]
        r = dist[inside]
        qq = np.broadcast_to(charges[i] * charges[:, None], inside.shape)[inside]
        kernel = r ** -3 - np.sum(weights * np.exp(-(r[:, None] / widths) ** 2), axis=1)
        pressure += np.einsum("p,pa,pb->ab", qq * kernel, d, d)
    return pressure / (2.0 * np.linalg.det(cell))

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
        compute_short_range_pressure({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_short_range_pressure({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark triclinic snapshot; the cutoff exceeds half the cell widths.
        {
            "setup": benchmark,
            "call": "compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0)",
            "gold_call": "_oracle_compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0)",
        },
        # Boundary: a single ion in a cubic cell shorter than the cutoff, so the whole
        # tensor comes from the ion's own periodic images.
        {
            "setup": "import numpy as np\ncell = 7.0 * np.eye(3)\npositions = np.array([[1.0, 2.0, 3.0]])\ncharges = np.array([2.0])\n",
            "call": "compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0)",
            "gold_call": "_oracle_compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 3.0, 10, 9.0)",
        },
        # Edge: ion pair with one ion outside the reference cell and a short cutoff in a
        # sheared cell, so only one image pair contributes.
        {
            "setup": "import numpy as np\ncell = np.array([[12.0, 3.0, 0.0], [0.0, 11.0, -2.0], [0.0, 0.0, 13.0]])\npositions = np.array([[0.5, 0.5, 0.5], [14.0, 2.5, -0.8]])\ncharges = np.array([1.0, -1.0])\n",
            "call": "compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 1.5, 8, 4.5)",
            "gold_call": "_oracle_compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 1.6, 1.5, 8, 4.5)",
        },
        # Edge: benchmark snapshot with a coarse splitting and fewer retained terms.
        {
            "setup": benchmark,
            "call": "compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 2.0, 4.0, 4, 9.0)",
            "gold_call": "_oracle_compute_short_range_pressure(positions.copy(), charges.copy(), cell.copy(), 2.0, 4.0, 4, 9.0)",
        },
        # Invalid: left-handed cell (negative determinant).
        invalid_case("np.zeros((2, 3)), np.array([1.0, -1.0]), np.diag([10.0, 10.0, -10.0]), 1.6, 3.0, 10, 9.0"),
        # Invalid: charges do not match the number of positions.
        invalid_case("np.zeros((2, 3)), np.array([1.0]), 10.0 * np.eye(3), 1.6, 3.0, 10, 9.0"),
    ]
