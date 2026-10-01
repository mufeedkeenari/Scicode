"""Whole-pipeline integration tests for the random-batch SOG shear-pressure task."""


def test_cases():
    """Return list of integration test case specifications."""
    benchmark = """import numpy as np
cell = np.array([[9.6, 1.4, -0.7], [0.0, 9.1, 1.1], [0.0, 0.0, 10.2]])
frac = np.array([[0.02, 0.03, 0.01], [0.47, 0.55, 0.04], [0.53, 0.06, 0.46], [0.05, 0.49, 0.52],
                 [0.51, 0.02, 0.07], [0.03, 0.54, 0.03], [0.06, 0.03, 0.55], [0.48, 0.46, 0.51]])
positions = frac @ cell.T
charges = np.array([1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0])
"""
    return [
        # Benchmark of the prompt: relative standard error of P_xy at batch size 128.
        {
            "setup": benchmark,
            "call": "compute_pressure_relative_error(positions, charges, cell, 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)",
            "gold_call": "0.37851954139827876",
            "tol": 1e-8,
            "final": True,
        },
        # Charge conjugation: reversing every charge leaves the pressure tensor (quadratic
        # in the charges) and the estimator variance (quartic) unchanged, so the relative
        # error must equal the benchmark value.
        {
            "setup": benchmark,
            "call": "compute_pressure_relative_error(positions, -charges, cell, 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)",
            "gold_call": "0.37851954139827876",
            "tol": 1e-8,
        },
        # Rigid translation of every ion by an arbitrary vector leaves all interparticle
        # separations and every |rho(k)|^2 unchanged, so the result must equal the benchmark.
        {
            "setup": benchmark + "shifted = positions + np.array([1.3, -2.2, 0.7])\n",
            "call": "compute_pressure_relative_error(shifted, charges, cell, 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)",
            "gold_call": "0.37851954139827876",
            "tol": 1e-8,
        },
        # Independent 2:1 salt snapshot (six ions) in a different triclinic cell with a
        # coarse splitting (b = 2, sigma = 4 A, six Gaussians, r_cut = 8 A), batch size 64,
        # shear component (y, z); expected value from the reference pipeline.
        {
            "setup": """import numpy as np
cell = np.array([[8.4, -1.1, 0.9], [0.0, 8.9, 1.6], [0.0, 0.0, 9.3]])
frac = np.array([[0.10, 0.12, 0.08], [0.62, 0.58, 0.55], [0.35, 0.70, 0.20],
                 [0.80, 0.30, 0.65], [0.15, 0.45, 0.90], [0.55, 0.05, 0.40]])
positions = frac @ cell.T
charges = np.array([2.0, -1.0, -1.0, 2.0, -1.0, -1.0])
""",
            "call": "compute_pressure_relative_error(positions, charges, cell, 2.0, 4.0, 6, 8.0, 6, 64, 1, 2)",
            "gold_call": "0.06056888896760986",
            "tol": 1e-8,
        },
    ]
