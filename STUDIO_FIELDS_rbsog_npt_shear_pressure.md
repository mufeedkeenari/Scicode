# Studio fields: Material Science -- Molecular Modeling (rbsog_npt_shear_pressure)

## 1 · Task, subject and area

**Task Name**

```text
Material Science -- Molecular Modeling
```

**Subject**

```text
Material Science
```

**Subject area / subdomain**

```text
Molecular Modeling
```

## 2 · Paper

**Source paper DOI**

```text
10.48550/arXiv.2602.23582
```

**Public open full-text URL (version used: v1)**

```text
https://arxiv.org/abs/2602.23582
```

**Source paper PDF to upload**

```text
arXiv 2602.23582v1 (Random batch sum-of-Gaussians method for molecular dynamics simulation of particle systems in the NPT ensemble), posted 27 Feb 2026
```

## 3 · Prompt and background

**Problem statement (the task)**

```markdown
Constant-pressure (NPT) molecular dynamics of electrolytes, ionic liquids and membranes needs the instantaneous Coulomb pressure tensor at every step, yet FFT-based Ewald solvers make its long-range part communication-bound at scale and their real-space cutoff discontinuities distort volume fluctuations. A recent stochastic scheme for constant-pressure simulations splits the pressure-related $1/r^3$ kernel with a bilateral-series sum of Gaussians into a short-range part truncated at a cutoff and a smooth long-range part evaluated in Fourier space, and it replaces the full Fourier sum by importance-sampled mini-batches of modes, using separate proposals for the radial and non-radial parts of the pressure tensor, with the non-radial modes obtained by a Metropolis–Hastings re-weighting of the radial samples. Given the charges, positions and cell tensor, the method returns the split pressure tensor and an unbiased stochastic estimate of it whose variance sets the batch size needed for stable dynamics.

Consider a snapshot of eight monovalent ions in a periodic triclinic cell whose matrix $h$ has the lattice vectors as its columns, $h_1 = (9.6, 0, 0)$ Å, $h_2 = (1.4, 9.1, 0)$ Å and $h_3 = (-0.7, 1.1, 10.2)$ Å, with Cartesian positions $r = h\,s$ for fractional coordinates $s$. The cations ($q = +1\,e$) sit at $s = (0.02, 0.03, 0.01)$, $(0.47, 0.55, 0.04)$, $(0.53, 0.06, 0.46)$ and $(0.05, 0.49, 0.52)$, and the anions ($q = -1\,e$) at $(0.51, 0.02, 0.07)$, $(0.03, 0.54, 0.03)$, $(0.06, 0.03, 0.55)$ and $(0.48, 0.46, 0.51)$. Use Gaussian units with a unit Coulomb prefactor (pressure in $e^2/\text{Å}^4$) and tinfoil boundary conditions. Split the $1/r^3$ kernel with the method's pressure-kernel parameters $\tilde b = 1.6$, $\tilde\sigma = 3.0$ Å and $\tilde M = 10$ retained long-range Gaussians and a real-space cutoff $r_c = 9.0$ Å, with the weight of the narrowest retained Gaussian rescaled so that the short-range pressure kernel built from these $\tilde M$ truncated terms is $C^0$-continuous at $r_c$; include every periodic image within $r_c$ in the real-space sum and converge all Fourier sums.

Take the deterministic sum-of-Gaussians value of the shear component $P_{xy}$ (short-range part plus full long-range part) as the reference, and consider the method's random-batch estimate of $P_{xy}$ with batch size $P = 128$, treating the 128 non-radial modes as independent draws from the method's non-radial importance distribution (an ideally mixed re-weighting chain). Report the rescaling factor of the narrowest Gaussian weight, the short-range and long-range contributions to $P_{xy}$, and the single-mode variance of the non-radial random-batch estimate of $P_{xy}$. The final answer is the relative standard error of the $P = 128$ estimate of $P_{xy}$, that is, its standard deviation divided by $|P_{xy}|$.
```

**Scientific background**

```markdown
Electrostatic interactions govern the structure, transport and thermodynamics of molten salts, ionic liquids, concentrated electrolytes and charged soft materials such as lipid membranes. In constant-pressure molecular dynamics, the barostat is driven at every time step by the instantaneous pressure (virial) tensor, so its Coulomb contribution must be accurate, smooth and cheap. Discontinuities introduced when a real-space interaction is truncated at a cutoff produce artificial pressure jumps and distorted cell-volume fluctuations, and under anisotropic or semi-isotropic pressure coupling they can corrupt membrane area and thickness fluctuations. Off-diagonal (shear) components of the pressure tensor are especially delicate, because they are small differences of large contributions and enter both anisotropic barostats and Green–Kubo viscosity calculations.

Periodic Coulomb lattice sums converge only conditionally and are traditionally evaluated with Ewald-type splittings, whose long-range part is computed in reciprocal space, usually with fast Fourier transforms. At large scale the global communication required by these transforms limits parallel efficiency. Alternative kernel decompositions represent long-range kernels as sums of Gaussians, which are separable and smooth, and stochastic random mini-batch methods replace full reciprocal-space sums by importance-sampled subsets of Fourier modes, giving linear cost and low communication. The usefulness of such stochastic schemes hinges on the variance of their unbiased estimates, which determines how many sampled modes are needed for stable and accurate dynamics.
```

## 4 · Browsing sources (Source 1)

**Query**

```text
random batch sum-of-Gaussians NPT ensemble pressure tensor 1/r^3 kernel splitting radial non-radial importance sampling measure recalibration
```

**Justification**

```markdown
The solver needs the method's pressure-specific constructions, which are not standard Ewald or energy-splitting results: the bilateral-series sum-of-Gaussians approximation of the pressure-related $1/r^3$ kernel with its weights and widths for the parameters $\tilde b,\tilde\sigma$ and its $\tilde M$-term truncation; the split of the long-range pressure tensor into a radial part proportional to the identity and a non-radial part along $k\otimes k$; and the non-radial importance distribution over Fourier modes, with its $|k|^4$ and $\tilde s_\ell^7$ weighting, that defines the single-mode random-batch estimate of the non-radial pressure. These choices fix the real-space kernel, the long-range tensor and the estimator variance; they cannot be inferred from generic Ewald or random-batch Ewald formulations, which use a different splitting and a single proposal.
```

**Target Source**

```text
https://arxiv.org/abs/2602.23582
```

## 5 · Golden solution – scientific reasoning chain

**Golden solution**

```markdown
<reasoning>

The calculation follows the random-batch sum-of-Gaussians (RBSOG) method for NPT simulations (arXiv:2602.23582): an SOG splitting of the pressure kernel $1/r^3$, the Fourier representation of its long-range part, and the radial/non-radial importance sampling of Fourier modes. Units: charges in $e$, lengths in Å, unit Coulomb prefactor, so pressures are in $e^2/\text{Å}^4$. The cell volume is $V=\det h = 891.072$ Å$^3$.

**Step 1 – SOG of the pressure kernel.** Specialising the bilateral series approximation of $r^{-\beta}$ (Eq. 2.16) to $\beta=3$ gives (Eq. 3.1)
$$\frac{1}{r^3}\approx\sum_{\ell}\tilde w_\ell\,e^{-r^2/\tilde s_\ell^2},\qquad \tilde w_\ell=\Big(\frac{\pi}{2}\Big)^{-1/2}\frac{\ln\tilde b}{\tilde b^{3\ell}\tilde\sigma^3},\qquad \tilde s_\ell=\sqrt2\,\tilde b^{\ell}\tilde\sigma .$$
The long-range kernel keeps $\ell=0,\dots,\tilde M-1$ (Eq. 3.3): $\tilde{\mathcal F}(r)=\sum_{\ell=0}^{9}\tilde w_\ell e^{-r^2/\tilde s_\ell^2}$. With $\tilde b=1.6$, $\tilde\sigma=3.0$ Å: $\tilde w_0=1.38892\times10^{-2}$ Å$^{-3}$, $\tilde s_0=4.24264$ Å, up to $\tilde s_9=291.552$ Å.

**Step 2 – $C^0$ rescaling of the narrowest Gaussian.** Eqs. (3.13)–(3.14) replace $\tilde w_0\to\tilde\omega\tilde w_0$ so that the short-range kernel $\tilde{\mathcal N}(r)=r^{-3}-\tilde{\mathcal F}(r)$ vanishes at $r_c$ for the truncated ten-term series:
$$\tilde\omega=\frac{r_c^{-3}-\sum_{\ell=1}^{9}\tilde w_\ell e^{-r_c^2/\tilde s_\ell^2}}{\tilde w_0\,e^{-r_c^2/\tilde s_0^2}}=1.011392 ,$$
so the rescaled narrowest weight is $1.404744\times10^{-2}$ Å$^{-3}$.

**Step 3 – Short-range pressure.** Inserting the splitting into the Coulomb pressure tensor (Eq. 2.14) gives the real-space part (Eq. 3.5)
$$P^{\mathcal N}=\frac{1}{2V}\sum_{n}{}'\sum_{i,j}q_iq_j\,\tilde{\mathcal N}(|r_{ij}+hn|)\,(r_{ij}+hn)\otimes(r_{ij}+hn),$$
with $\tilde{\mathcal N}=0$ beyond $r_c$. The perpendicular cell widths are 9.45, 9.05 and 10.20 Å, so $r_c=9$ Å exceeds half of every width: 190 ordered pair–image terms lie inside $r_c$, 134 of them beyond the minimum image, and all must be summed. The result is
$$P^{\mathcal N}_{xy}=-8.8851\times10^{-5}\ e^2/\text{Å}^4 .$$

**Step 4 – Long-range pressure in Fourier space.** With $\hat f(k)=\int f(r)e^{-ik\cdot r}dr$, the Fourier transform of $\tilde{\mathcal F}(|r|)\,r\otimes r$ is $\sum_\ell\tfrac12\pi^{3/2}\tilde w_\ell\tilde s_\ell^5e^{-\tilde s_\ell^2k^2/4}\big(I-\tfrac{\tilde s_\ell^2}{2}k\otimes k\big)$. Poisson summation over the lattice $k=2\pi h^{-\top}m$, dropping $k=0$ (tinfoil, neutral system), gives Theorem 1 (Eq. 3.8):
$$P^{\mathcal F}=\frac{\pi^{3/2}}{4V^2}\sum_{k\neq0}\sum_{\ell=0}^{9}\tilde w_\ell\tilde s_\ell^5e^{-\tilde s_\ell^2k^2/4}\Big(I-\frac{\tilde s_\ell^2}{2}k\otimes k\Big)|\rho(k)|^2,\qquad \rho(k)=\sum_jq_je^{2\pi i\,m\cdot s_j}.$$
This splits into the radial part $P^{\mathcal F,\mathrm r}$ (Eq. 3.17), which is proportional to $I$ ($6.0989\times10^{-5}\,I$), and the non-radial part $P^{\mathcal F,\mathrm{nr}}=-\frac{\pi^{3/2}}{8V^2}\sum_{k\ne0}\sum_\ell\tilde w_\ell\tilde s_\ell^7e^{-\tilde s_\ell^2k^2/4}|\rho|^2k\otimes k$ (Eq. 3.18). The sum is converged at $|m_d|\le4$, since the narrowest Gaussian factor $e^{-\tilde s_0^2k^2/4}$ is below $3\times10^{-19}$ for every mode outside that set. Then
$$P^{\mathcal F}_{xy}=P^{\mathcal F,\mathrm{nr}}_{xy}=-7.6565\times10^{-5}\ e^2/\text{Å}^4,\qquad P_{xy}=P^{\mathcal N}_{xy}+P^{\mathcal F}_{xy}=-1.65417\times10^{-4}\ e^2/\text{Å}^4 .$$
This SOG value agrees with an independent Ewald evaluation of the exact Coulomb $P_{xy}$ ($-1.6587\times10^{-4}$) to 0.3 %, consistent with the decomposition accuracy of $\tilde b=1.6$.

**Step 5 – Which estimator carries the shear noise.** The radial random-batch estimate (Eq. 3.21) is $P^{\mathcal F,\mathrm r*}=\frac{S^{\mathrm r}}{4V^2P}\sum_p|\rho(k_p)|^2/|k_p|^2\,I$, which is isotropic for any sample, so it adds no noise to $P_{xy}$. The shear noise comes entirely from the non-radial estimate. Modes are drawn from the non-radial proposal (Eqs. 3.22–3.23)
$$\mathscr P^{\mathrm{nr}}(k)=\frac{1}{S^{\mathrm{nr}}}\sum_\ell\pi^{3/2}\tilde w_\ell\tilde s_\ell^7|k|^4e^{-\tilde s_\ell^2k^2/4},\qquad S^{\mathrm{nr}}=892.989,$$
which measure recalibration reaches as the stationary target of the chain built from the radial samples. Each mode gives the unbiased single-mode estimate (Eq. 3.24 with $P=1$)
$$X(k)=-\frac{S^{\mathrm{nr}}|\rho(k)|^2k_xk_y}{8V^2|k|^4},\qquad \mathbb E[X]=P^{\mathcal F,\mathrm{nr}}_{xy}.$$

**Step 6 – Exact single-mode variance.** For an importance-sampled lattice sum $\mu=\sum_kf(k)$ estimated by $f(k)/\mathscr P(k)$, the variance is $\sum_kf(k)^2/\mathscr P(k)-\mu^2$. Here this gives
$$\mathrm{Var}[X]=\frac{\pi^{3/2}S^{\mathrm{nr}}}{64V^4}\sum_{k\ne0}\Big(\sum_\ell\tilde w_\ell\tilde s_\ell^7e^{-\tilde s_\ell^2k^2/4}\Big)\frac{|\rho(k)|^4k_x^2k_y^2}{|k|^4}-\big(P^{\mathcal F,\mathrm{nr}}_{xy}\big)^2=5.0182\times10^{-7}\ e^4/\text{Å}^8 .$$
Proposition 7 (Eq. 3.43) is an inequality whose non-radial prefactor $1/36$ exceeds the exact $1/64$; used as an equality, it would overstate the variance at $8.97\times10^{-7}$.

**Step 7 – Batch of 128 modes.** The batch estimate averages 128 independent single-mode estimates, so its standard deviation is
$$\sqrt{\mathrm{Var}[X]/128}=6.2614\times10^{-5}\ e^2/\text{Å}^4 .$$

**Step 8 – Relative standard error.**
$$\frac{6.2614\times10^{-5}}{|{-1.65417\times10^{-4}}|}=0.37852 .$$
At $P=128$ the random-batch noise in the instantaneous shear pressure of this small eight-ion cell is therefore about 38 % of its deterministic value.

</reasoning>

<final_answer>0.3785195413982788</final_answer>
```

## 6 · Rubric (11 items)

### Rubric item 1

**Item 1 · Criterion**

```markdown
Identifies the sum-of-Gaussians approximation of the pressure kernel $1/r^3$, with splitting parameters $\tilde b,\tilde\sigma$ and $\tilde M$ retained long-range Gaussians, as $\sum_\ell\tilde w_\ell e^{-r^2/\tilde s_\ell^2}$ with $\tilde w_\ell=(\pi/2)^{-1/2}\ln\tilde b/(\tilde b^{3\ell}\tilde\sigma^3)$ and $\tilde s_\ell=\sqrt2\,\tilde b^{\ell}\tilde\sigma$, its long-range part keeping $\ell=0,\dots,\tilde M-1$.
```

**Item 1 · Category**

```text
Browsing
```

**Item 1 · Weight**

```text
3
```

**Item 1 · Description**

```markdown
This is the bilateral-series approximation of $r^{-\beta}$ at $\beta=3$, with the narrowest retained width $\sqrt2\tilde\sigma=4.243$ Å for $\tilde\sigma=3$ Å. Algebraically equivalent prefactors, such as $\sqrt2\ln\tilde b/(\sqrt\pi\,\tilde\sigma^3)\,\tilde b^{-3\ell}$, count.
```

**Item 1 · Source**

```text
https://arxiv.org/abs/2602.23582
```

### Rubric item 2

**Item 2 · Criterion**

```markdown
Identifies the non-radial importance distribution over nonzero Fourier wavevectors $k$, with $\tilde w_\ell,\tilde s_\ell$ the weights and widths of the retained long-range Gaussians of the $1/r^3$ splitting ($\ell=0,\dots,\tilde M-1$), as $\mathscr P^{\rm nr}(k)\propto\sum_\ell\tilde w_\ell\tilde s_\ell^{7}|k|^4e^{-\tilde s_\ell^2|k|^2/4}$.
```

**Item 2 · Category**

```text
Browsing
```

**Item 2 · Weight**

```text
5
```

**Item 2 · Description**

```markdown
The method samples the non-radial pressure part from this proposal, which carries an extra factor $\tilde s_\ell^2|k|^2$ relative to the radial proposal $\propto\sum_\ell\tilde w_\ell\tilde s_\ell^5|k|^2e^{-\tilde s_\ell^2|k|^2/4}$. Any normalization constant is acceptable. A proposal proportional to the magnitude of the $P_{xy}$ summand, or the radial proposal, does not count.
```

**Item 2 · Source**

```text
https://arxiv.org/abs/2602.23582
```

### Rubric item 3

**Item 3 · Criterion**

```markdown
Identifies the radial random-batch estimate of the long-range pressure tensor as an isotropic tensor (a multiple of the identity) that adds zero variance to the estimate of $P_{xy}$.
```

**Item 3 · Category**

```text
Browsing
```

**Item 3 · Weight**

```text
2
```

**Item 3 · Description**

```markdown
The method's radial part collects the identity-proportional term of each Fourier mode, so each of its samples is isotropic. The shear-component noise therefore comes only from the non-radial estimate.
```

**Item 3 · Source**

```text
https://arxiv.org/abs/2602.23582
```

### Rubric item 4

**Item 4 · Criterion**

```markdown
Derives the long-range (Fourier-space) pressure tensor, with $V=\det h$, $\tilde w_\ell,\tilde s_\ell$ the weights and widths of the retained long-range Gaussians of the $1/r^3$ splitting ($\ell=0,\dots,\tilde M-1$), wavevectors $k\neq0$ and charge structure factor $\rho(k)=\sum_jq_je^{ik\cdot r_j}$, as $P^{\mathcal F}=\frac{\pi^{3/2}}{4V^2}\sum_{k\ne0}\sum_\ell\tilde w_\ell\tilde s_\ell^5e^{-\tilde s_\ell^2|k|^2/4}\big(I-\tfrac{\tilde s_\ell^2}{2}k\otimes k\big)|\rho(k)|^2$.
```

**Item 4 · Category**

```text
Scientific Reasoning
```

**Item 4 · Weight**

```text
4
```

**Item 4 · Description**

```markdown
This follows from Poisson summation of the Fourier transform of the long-range kernel times $r\otimes r$, dropping $k=0$ for tinfoil boundary conditions and a neutral system. Equivalent forms, such as writing the radial and non-radial parts separately or summing over integer vectors $m$, count.
```

### Rubric item 5

**Item 5 · Criterion**

```markdown
Identifies the Fourier wavevectors of the triclinic cell, with $h$ the matrix whose columns are the lattice vectors and $m$ integer vectors, as $k=2\pi(h^{-1})^{\top}m$.
```

**Item 5 · Category**

```text
Scientific Reasoning
```

**Item 5 · Weight**

```text
2
```

**Item 5 · Description**

```markdown
This makes $k\cdot r=2\pi\,m\cdot s$ for fractional coordinates $s$ and reproduces the periodicity of the sheared cell. Equivalent statements count, such as $k=\sum_i m_ib_i$ with reciprocal vectors $a_i\cdot b_j=2\pi\delta_{ij}$, or $k=2\pi h^{-1}m$ when $h$ is defined with the lattice vectors as rows. With the lattice vectors as the columns of $h$, $2\pi h^{-1}m$ gives wrong wavevectors for this non-orthogonal cell.
```

### Rubric item 6

**Item 6 · Criterion**

```markdown
Computes the rescaling factor $\tilde\omega$ of the narrowest retained Gaussian weight as 1.0114 (accepting 1.0112 to 1.0116).
```

**Item 6 · Category**

```text
Scientific Reasoning
```

**Item 6 · Weight**

```text
4
```

**Item 6 · Description**

```markdown
Imposing $1/r_c^3=\sum_{\ell=0}^{9}\tilde w_\ell e^{-r_c^2/\tilde s_\ell^2}$ at $r_c=9$ Å, with $\tilde w_0\to\tilde\omega\tilde w_0$, gives 1.011392. The range also accepts 1.011303, obtained with an untruncated tail, which changes the final answer by less than 0.01 %.
```

### Rubric item 7

**Item 7 · Criterion**

```markdown
Computes the short-range (real-space) contribution to $P_{xy}$, summed over all periodic images within $r_c=9$ Å, as $-8.885\times10^{-5}\ e^2/\text{Å}^4$ (within ±0.5 %).
```

**Item 7 · Category**

```text
Scientific Reasoning
```

**Item 7 · Weight**

```text
5
```

**Item 7 · Description**

```markdown
The rescaled short-range kernel $1/r^3-\sum_\ell\tilde w_\ell e^{-r^2/\tilde s_\ell^2}$ is summed over pair and image separations below $r_c$. Because $r_c$ exceeds half of every perpendicular cell width (9.45, 9.05, 10.20 Å), a minimum-image sum ($\approx-4.79\times10^{-5}$) is wrong. The same value expressed in other pressure units counts.
```

### Rubric item 8

**Item 8 · Criterion**

```markdown
Computes the long-range (Fourier-space) contribution to $P_{xy}$ as $-7.657\times10^{-5}\ e^2/\text{Å}^4$ (within ±0.5 %).
```

**Item 8 · Category**

```text
Scientific Reasoning
```

**Item 8 · Weight**

```text
5
```

**Item 8 · Description**

```markdown
This is the off-diagonal entry of the converged Fourier-space pressure tensor; integer vectors with $|m_d|\le4$ already converge it, and only the non-radial part contributes. Together with the short-range part it gives $P_{xy}=-1.6542\times10^{-4}\ e^2/\text{Å}^4$. The same value expressed in other pressure units counts.
```

### Rubric item 9

**Item 9 · Criterion**

```markdown
States the single-sample variance of an importance-sampled sum $\mu=\sum_kf(k)$, estimated by $f(k)/\mathscr P(k)$ with $k$ drawn from $\mathscr P$, as $\sum_kf(k)^2/\mathscr P(k)-\mu^2$.
```

**Item 9 · Category**

```text
Scientific Reasoning
```

**Item 9 · Weight**

```text
3
```

**Item 9 · Description**

```markdown
This exact identity, applied with $f$ the $P_{xy}$ summand of the non-radial tensor and $\mathscr P=\mathscr P^{\rm nr}$, gives the requested variance. An equivalent form that makes the sampling distribution explicit, such as $\mathbb E_{k\sim\mathscr P}\big[(f(k)/\mathscr P(k))^2\big]-\mu^2$ for the non-radial estimate, counts.
```

### Rubric item 10

**Item 10 · Criterion**

```markdown
Computes the single-mode variance of the non-radial random-batch estimate of $P_{xy}$ as $5.018\times10^{-7}\ e^4/\text{Å}^8$ (within ±1 %).
```

**Item 10 · Category**

```text
Scientific Reasoning
```

**Item 10 · Weight**

```text
7
```

**Item 10 · Description**

```markdown
The single-mode estimate is $-S^{\rm nr}|\rho(k)|^2k_xk_y/(8V^2|k|^4)$ with mean $-7.657\times10^{-5}$, where $S^{\rm nr}=\pi^{3/2}\sum_{k\ne0}\sum_\ell\tilde w_\ell\tilde s_\ell^7|k|^4e^{-\tilde s_\ell^2|k|^2/4}$ normalizes $\mathscr P^{\rm nr}$ ($\tilde w_\ell,\tilde s_\ell$ the long-range Gaussian weights, with $\tilde w_0$ rescaled by $\tilde\omega$, and widths). A value of $8.97\times10^{-7}$, from treating the source's variance bound as an equality, or $9.09\times10^{-7}$, from sampling the radial proposal, does not count. The same value expressed in other units of squared pressure counts.
```

### Rubric item 11

**Item 11 · Criterion**

```markdown
Computes the relative standard error of the random-batch estimate of $P_{xy}$ with batch size $P=128$ as 0.3785 (within ±0.5 %, i.e. 0.3766 to 0.3804).
```

**Item 11 · Category**

```text
Scientific Reasoning
```

**Item 11 · Weight**

```text
8
```

**Item 11 · Description**

```markdown
This equals $\sqrt{5.018\times10^{-7}/128}/|{-1.6542\times10^{-4}}|=6.261\times10^{-5}/1.6542\times10^{-4}$. Omitting the narrowest-weight rescaling gives 0.3748, and treating the variance bound as an equality gives 0.506; neither counts.
```

## 7 · Sub-problems (8 steps; step 08 is the final orchestrator)

### Step 01 · compute_pressure_kernel_gaussians

**Step 01 · Step name**

```text
compute_pressure_kernel_gaussians
```

**Step 01 · Description**

```text
Compute the weights and widths of the retained long-range Gaussian terms in the bilateral-series sum-of-Gaussians approximation of the pressure-related 1/r^3 kernel.
```

**Step 01 · Scientific background**

```text
The instantaneous Coulomb pressure tensor of a periodic charge system is governed by the
radial kernel 1/r^3. A smooth short-/long-range splitting of this kernel represents its
long-range part as a finite sum of Gaussians, each written as w_l * exp(-r^2 / s_l^2) with
weight w_l and width s_l (lengths in Angstrom, weights in Angstrom^-3). The splitting is
controlled by the bilateral-series parameters b > 1 and sigma > 0 (Angstrom) and by the
number n_terms of retained long-range Gaussians.

Return the retained terms ordered from the narrowest (index 0) to the widest
(index n_terms - 1), before any rescaling of the narrowest weight.
```

**Step 01 · Signature (shown to the LLM)**

```python
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
```

**Step 01 · Expected return line**

```text
tuple (weights, widths) of two float np.ndarray of shape (n_terms,): the weights w_l in Angstrom^-3 (before rescaling of the narrowest weight) and the widths s_l in Angstrom of the terms w_l * exp(-r^2 / s_l^2), both ordered from the narrowest to the widest Gaussian
```

**Step 01 · Oracle (gold solution)**

```python
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
```

**Step 01 · Test cases**

```python
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
```

### Step 02 · compute_narrowest_weight_factor

**Step 02 · Step name**

```text
compute_narrowest_weight_factor
```

**Step 02 · Description**

```text
Compute the factor multiplying the weight of the narrowest retained Gaussian of the 1/r^3 kernel splitting so that the short-range pressure kernel is continuous at the real-space cutoff.
```

**Step 02 · Scientific background**

```text
The short-range part of the pressure kernel is truncated at a real-space cutoff r_cut. A
cutoff discontinuity in the instantaneous pressure produces artificial pressure jumps and
volume fluctuations in constant-pressure simulations, so the splitting is adjusted by
multiplying the weight of its narrowest retained Gaussian (index 0 of the retained terms)
by a dimensionless factor.

The continuity condition is imposed on the truncated long-range kernel built from exactly
the n_terms retained Gaussians, not on an untruncated series. The factor must be positive,
so that every long-range Gaussian keeps a positive weight.
```

**Step 02 · Signature (shown to the LLM)**

```python
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
```

**Step 02 · Expected return line**

```text
float, the dimensionless positive factor that multiplies the weight of the narrowest retained Gaussian
```

**Step 02 · Oracle (gold solution)**

```python
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
```

**Step 02 · Test cases**

```python
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
```

### Step 03 · compute_short_range_pressure

**Step 03 · Step name**

```text
compute_short_range_pressure
```

**Step 03 · Description**

```text
Compute the short-range (real-space) contribution to the instantaneous Coulomb pressure tensor of a periodic point-charge system under the rescaled sum-of-Gaussians splitting of the 1/r^3 kernel.
```

**Step 03 · Scientific background**

```text
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
```

**Step 03 · Signature (shown to the LLM)**

```python
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
```

**Step 03 · Expected return line**

```text
np.ndarray of shape (3, 3), the symmetric short-range Coulomb pressure tensor in e^2 / Angstrom^4, indexed (x, y, z) = (0, 1, 2)
```

**Step 03 · Oracle (gold solution)**

```python
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
```

**Step 03 · Test cases**

```python
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
```

### Step 04 · compute_structure_factor_power

**Step 04 · Step name**

```text
compute_structure_factor_power
```

**Step 04 · Description**

```text
Compute the squared modulus of the charge structure factor of a periodic point-charge system at a given list of reciprocal-lattice modes.
```

**Step 04 · Scientific background**

```text
The triclinic cell matrix has the lattice vectors as its columns, so a particle with
fractional coordinates s sits at r = cell @ s. Fourier modes of the periodic cell are
labelled by integer vectors m; mode m corresponds to the wavevector
k = 2 * pi * inv(cell).T @ m (Angstrom^-1), so that k . (cell @ n) is a multiple of 2 * pi
for every integer lattice translation n. The charge structure factor is the Fourier sum of
the charges over the particle positions, and its squared modulus (units of e^2) is
independent of the sign convention of the Fourier exponent and invariant under lattice
translations of individual particles.
```

**Step 04 · Signature (shown to the LLM)**

```python
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
```

**Step 04 · Expected return line**

```text
np.ndarray of shape (K,), the float squared modulus of the charge structure factor in e^2 at each mode, in the order of the rows of modes
```

**Step 04 · Oracle (gold solution)**

```python
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
```

**Step 04 · Test cases**

```python
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
```

### Step 05 · compute_long_range_pressure

**Step 05 · Step name**

```text
compute_long_range_pressure
```

**Step 05 · Description**

```text
Compute the long-range (Fourier-space) contribution to the instantaneous Coulomb pressure tensor under the rescaled sum-of-Gaussians splitting of the 1/r^3 kernel, separated into its radial and non-radial parts.
```

**Step 05 · Scientific background**

```text
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
```

**Step 05 · Signature (shown to the LLM)**

```python
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
```

**Step 05 · Expected return line**

```text
tuple (radial, nonradial) of two float np.ndarray of shape (3, 3) in e^2 / Angstrom^4: the radial part (a multiple of the identity tensor) and the symmetric non-radial part of the long-range pressure tensor
```

**Step 05 · Oracle (gold solution)**

```python
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
```

**Step 05 · Test cases**

```python
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
```

### Step 06 · compute_nonradial_normalization

**Step 06 · Step name**

```text
compute_nonradial_normalization
```

**Step 06 · Description**

```text
Compute the normalization constant of the non-radial importance-sampling distribution over Fourier modes used by the random-batch estimator of the non-radial long-range pressure.
```

**Step 06 · Scientific background**

```text
The random-batch method replaces the full Fourier sum of the non-radial long-range pressure
by an average over a small batch of modes drawn from a discrete importance distribution
defined by the rescaled long-range Gaussians. The distribution is supported on the
nonzero integer modes m with max(|m_x|, |m_y|, |m_z|) <= m_max, where mode m labels the
wavevector k = 2 * pi * inv(cell).T @ m for a cell matrix whose columns are the lattice
vectors. The returned constant is the normalization of that distribution, summed over
exactly this truncated mode set. The importance distribution depends only on the cell and
the splitting, not on the particle configuration.
```

**Step 06 · Signature (shown to the LLM)**

```python
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
```

**Step 06 · Expected return line**

```text
float, the dimensionless normalization constant of the non-radial Fourier-mode importance distribution over the truncated mode set
```

**Step 06 · Oracle (gold solution)**

```python
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
```

**Step 06 · Test cases**

```python
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
```

### Step 07 · compute_nonradial_variance

**Step 07 · Step name**

```text
compute_nonradial_variance
```

**Step 07 · Description**

```text
Compute the single-mode variance of the random-batch estimate of one Cartesian component of the non-radial long-range Coulomb pressure tensor when the Fourier mode is drawn from the non-radial importance distribution.
```

**Step 07 · Scientific background**

```text
A random-batch estimate of the non-radial long-range pressure averages an unbiased
single-mode estimate over P independently drawn modes, so its variance is the single-mode
variance divided by P. This step returns that single-mode (P = 1) variance for the (mu, nu)
component, with the mode drawn exactly from the non-radial importance distribution on the
nonzero integer modes with max(|m_x|, |m_y|, |m_z|) <= m_max (mode m labels the wavevector
k = 2 * pi * inv(cell).T @ m). The estimate's mean is the corresponding component of the
non-radial long-range tensor summed over the same mode set.

Units: charges in e, lengths in Angstrom, unit Coulomb prefactor; the variance is in
e^4 / Angstrom^8. Cartesian indices are (x, y, z) = (0, 1, 2).
```

**Step 07 · Signature (shown to the LLM)**

```python
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
```

**Step 07 · Expected return line**

```text
float, the single-mode variance of the (mu, nu) component of the non-radial random-batch estimate in e^4 / Angstrom^8
```

**Step 07 · Oracle (gold solution)**

```python
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
```

**Step 07 · Test cases**

```python
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
```

### Step 08 · compute_pressure_relative_error (FINAL / orchestrator)

**Step 08 · Step name**

```text
compute_pressure_relative_error
```

**Step 08 · Description**

```text
Compute the relative standard error of the random-batch estimate, under the sum-of-Gaussians splitting, of an off-diagonal (shear) component of the instantaneous Coulomb pressure tensor for a given batch size (end-to-end pipeline; orchestrator).
```

**Step 08 · Scientific background**

```text
The deterministic reference value of the (mu, nu) shear component is the rescaled
sum-of-Gaussians pressure tensor: its short-range real-space part plus its full long-range
Fourier-space part (radial plus non-radial) on the truncated mode set. In the random-batch
method the short-range part is evaluated exactly, while the long-range part is replaced by
importance-sampled mini-batches of Fourier modes. Here the batch_size modes of the
non-radial estimate are treated as independent draws from the non-radial importance
distribution, i.e. an ideally mixed Metropolis-Hastings re-weighting chain.

The returned quantity is the standard deviation of the batch estimate of the (mu, nu)
component divided by the magnitude of its deterministic reference value (dimensionless).
Only off-diagonal components (mu != nu) are accepted.
```

**Step 08 · Signature (shown to the LLM)**

```python
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
```

**Step 08 · Expected return line**

```text
float, the dimensionless relative standard error: the standard deviation of the batch estimate of the (mu, nu) pressure component divided by the magnitude of its deterministic sum-of-Gaussians value
```

**Step 08 · Oracle (gold solution)**

```python
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
```

**Step 08 · Test cases**

```python
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
```

## 7b · Integration tests – whole pipeline

**Integration tests**

```python
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
```

## 8 · Classification

**Difficulty**

```text
difficult
```

**Problem type**

```text
Numerical computation (multi-step scientific coding)
```

**Answer type**

```text
numeric
```
