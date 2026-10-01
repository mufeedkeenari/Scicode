"""
Compute the squared modulus of the charge structure factor of a periodic point-charge system at a given list of reciprocal-lattice modes.

The triclinic cell matrix has the lattice vectors as its columns, so a particle with
fractional coordinates s sits at r = cell @ s. Fourier modes of the periodic cell are
labelled by integer vectors m; mode m corresponds to the wavevector
k = 2 * pi * inv(cell).T @ m (Angstrom^-1), so that k . (cell @ n) is a multiple of 2 * pi
for every integer lattice translation n. The charge structure factor is the Fourier sum of
the charges over the particle positions, and its squared modulus (units of e^2) is
independent of the sign convention of the Fourier exponent and invariant under lattice
translations of individual particles.
"""

# =============================================================================
# FUNCTION SIGNATURE
# =============================================================================

def compute_structure_factor_power(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", modes: "np.ndarray") -> "np.ndarray":
    '''Squared modulus of the charge structure factor at integer reciprocal-lattice modes.

    Parameters
    ----------
    positions : np.ndarray
        Cartesian positions in Angstrom, shape (N, 3) with N >= 1. Not modified.
    charges : np.ndarray
        Charges in units of e, shape (N,). Not modified.
    cell : np.ndarray
        Cell matrix of shape (3, 3) whose columns are the lattice vectors in Angstrom;
        its determinant must be positive. Not modified.
    modes : np.ndarray
        Integer-valued mode vectors m, shape (K, 3) with K >= 1. Not modified.

    Returns
    -------
    power : np.ndarray
        Float array of shape (K,), the squared modulus of the charge structure factor in
        e^2 at each mode, in the order of the rows of modes.

    Raises
    ------
    ValueError
        If positions, charges or cell have invalid shapes or non-finite entries, if the
        cell determinant is not positive, or if modes is not a (K, 3) array of integer
        values with K >= 1.
    '''
    return power

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


def _oracle_compute_structure_factor_power(positions: "np.ndarray", charges: "np.ndarray", cell: "np.ndarray", modes: "np.ndarray") -> "np.ndarray":
    positions, charges, cell = _validate_system(positions, charges, cell)
    modes = np.array(modes, dtype=float)
    if modes.ndim != 2 or modes.shape[1] != 3 or modes.shape[0] < 1:
        raise ValueError("modes must have shape (K, 3) with K >= 1")
    if not np.all(np.isfinite(modes)) or np.any(modes != np.round(modes)):
        raise ValueError("modes must contain integer values")
    frac = positions @ np.linalg.inv(cell).T
    phase = 2.0 * np.pi * (modes @ frac.T)
    rho = np.exp(1j * phase) @ charges
    return rho.real ** 2 + rho.imag ** 2

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
r = np.arange(-2, 3)
modes = np.stack(np.meshgrid(r, r, r, indexing="ij"), axis=-1).reshape(-1, 3)
"""
    def invalid_case(args):
        setup = f"""import numpy as np
def run_model():
    try:
        compute_structure_factor_power({args})
        return 0
    except ValueError:
        return 1
def run_oracle():
    try:
        _oracle_compute_structure_factor_power({args})
        return 0
    except ValueError:
        return 1
"""
        return {"setup": setup, "call": "run_model()", "gold_call": "run_oracle()"}

    return [
        # Normal: benchmark snapshot on all modes with |m_d| <= 2 (including m = 0).
        {
            "setup": benchmark,
            "call": "compute_structure_factor_power(positions.copy(), charges.copy(), cell.copy(), modes.copy())",
            "gold_call": "_oracle_compute_structure_factor_power(positions.copy(), charges.copy(), cell.copy(), modes.copy())",
        },
        # Boundary: one ion and one mode; the power equals the squared charge.
        {
            "setup": "import numpy as np\ncell = np.diag([5.0, 6.0, 7.0])\npositions = np.array([[1.0, 2.0, 3.0]])\ncharges = np.array([-2.0])\nmodes = np.array([[1, -1, 2]])\n",
            "call": "compute_structure_factor_power(positions.copy(), charges.copy(), cell.copy(), modes.copy())",
            "gold_call": "_oracle_compute_structure_factor_power(positions.copy(), charges.copy(), cell.copy(), modes.copy())",
        },
        # Edge: strongly sheared cell with ions outside the reference cell; distinguishes
        # the inv(cell).T mode convention from inv(cell).
        {
            "setup": "import numpy as np\ncell = np.array([[8.0, 4.0, 2.0], [0.0, 7.0, 3.5], [0.0, 0.0, 9.0]])\npositions = np.array([[-3.0, 1.0, 12.0], [2.0, 9.5, -1.0], [5.0, -2.0, 4.0]])\ncharges = np.array([2.0, -1.0, -1.0])\nmodes = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, -2, 1], [-3, 1, 2]])\n",
            "call": "compute_structure_factor_power(positions.copy(), charges.copy(), cell.copy(), modes.copy())",
            "gold_call": "_oracle_compute_structure_factor_power(positions.copy(), charges.copy(), cell.copy(), modes.copy())",
        },
        # Invalid: non-integer mode vector.
        invalid_case("np.zeros((1, 3)), np.array([1.0]), 10.0 * np.eye(3), np.array([[0.5, 0.0, 1.0]])"),
        # Invalid: modes with the wrong shape.
        invalid_case("np.zeros((1, 3)), np.array([1.0]), 10.0 * np.eye(3), np.array([1, 0, 0])"),
    ]
