## Source 1
**Query:** random batch sum-of-Gaussians NPT ensemble pressure tensor 1/r^3 kernel splitting radial non-radial importance sampling measure recalibration

**Justification:** The solver needs the method's pressure-specific constructions for the long-range part, which are not standard Ewald or energy-splitting results: the split of the long-range pressure tensor into a radial part proportional to the identity and a non-radial part along $k\otimes k$, with the radial random-batch estimate isotropic; and the non-radial importance distribution over Fourier modes, with its $|k|^4$ and $\tilde s_\ell^7$ weighting, that defines the single-mode random-batch estimate of the non-radial pressure. These choices fix which part of the estimate carries the shear noise and the estimator variance; they cannot be inferred from generic Ewald or random-batch Ewald formulations, which use a different kernel splitting and a single proposal.

**Target Source:** https://arxiv.org/abs/2602.23582
