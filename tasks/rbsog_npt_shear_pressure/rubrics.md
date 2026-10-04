### 1

#### Criterion

Identifies the sum-of-Gaussians approximation of the pressure kernel $1/r^3$, with splitting parameters $\tilde b,\tilde\sigma$ and $\tilde M$ retained long-range Gaussians, as $\sum_\ell\tilde w_\ell e^{-r^2/\tilde s_\ell^2}$ with $\tilde w_\ell=(\pi/2)^{-1/2}\ln\tilde b/(\tilde b^{3\ell}\tilde\sigma^3)$ and $\tilde s_\ell=\sqrt2\,\tilde b^{\ell}\tilde\sigma$, its long-range part keeping $\ell=0,\dots,\tilde M-1$.

#### Category

Scientific Reasoning

#### Weight

5

#### Description

This is the bilateral-series approximation of $r^{-\beta}$ at $\beta=3$, with the narrowest retained width $\sqrt2\tilde\sigma=4.243$ Å for $\tilde\sigma=3$ Å. Algebraically equivalent forms count, such as the prefactor $\sqrt2\ln\tilde b/(\sqrt\pi\,\tilde\sigma^3)\,\tilde b^{-3\ell}$ or the standard-deviation form $\tilde w_\ell e^{-r^2/(2\sigma_\ell^2)}$ with $\sigma_\ell=\tilde b^{\ell}\tilde\sigma$ and the same weights. A series with a different prefactor (for example the $1/r$ weights $(\pi/2)^{-1/2}\ln\tilde b/(\tilde b^{\ell}\tilde\sigma)$), Gaussians $e^{-r^2/(\tilde b^{2\ell}\tilde\sigma^2)}$ without the factor 2 in the exponent, or a retained range other than $\ell=0,\dots,\tilde M-1$ does not count.

### 2

#### Criterion

Identifies the non-radial importance distribution over nonzero Fourier wavevectors $k$, with $\tilde w_\ell,\tilde s_\ell$ the weights and widths of the retained long-range Gaussians of the $1/r^3$ splitting ($\ell=0,\dots,\tilde M-1$), as $\mathscr P^{\rm nr}(k)\propto\sum_\ell\tilde w_\ell\tilde s_\ell^{7}|k|^4e^{-\tilde s_\ell^2|k|^2/4}$.

#### Category

Browsing

#### Weight

10

#### Description

The method samples the non-radial pressure part from this proposal, which carries an extra factor $\tilde s_\ell^2|k|^2$ relative to the radial proposal $\propto\sum_\ell\tilde w_\ell\tilde s_\ell^5|k|^2e^{-\tilde s_\ell^2|k|^2/4}$. Any normalization constant is acceptable. A proposal proportional to the magnitude of the $P_{xy}$ summand, or the radial proposal, does not count.

#### Source

https://arxiv.org/abs/2602.23582

### 3

#### Criterion

Identifies the radial random-batch estimate of the long-range pressure tensor as an isotropic tensor (a multiple of the identity) that adds zero variance to the estimate of $P_{xy}$.

#### Category

Browsing

#### Weight

2

#### Description

The method's radial part collects the identity-proportional term of each Fourier mode, so each of its samples is isotropic. The shear-component noise therefore comes only from the non-radial estimate. A response that attributes any part of the $P_{xy}$ noise to the radial estimate, or adds a radial variance term to the shear variance, does not count.

#### Source

https://arxiv.org/abs/2602.23582

### 4

#### Criterion

Derives the long-range (Fourier-space) pressure tensor, with $V=\det h$, $\tilde w_\ell,\tilde s_\ell$ the weights and widths of the retained long-range Gaussians of the $1/r^3$ splitting ($\ell=0,\dots,\tilde M-1$), wavevectors $k\neq0$ and charge structure factor $\rho(k)=\sum_jq_je^{ik\cdot r_j}$, as $P^{\mathcal F}=\frac{\pi^{3/2}}{4V^2}\sum_{k\ne0}\sum_\ell\tilde w_\ell\tilde s_\ell^5e^{-\tilde s_\ell^2|k|^2/4}\big(I-\tfrac{\tilde s_\ell^2}{2}k\otimes k\big)|\rho(k)|^2$.

#### Category

Scientific Reasoning

#### Weight

5

#### Description

This follows from Poisson summation of the Fourier transform of the long-range kernel times $r\otimes r$, dropping $k=0$ for tinfoil boundary conditions and a neutral system. Equivalent forms, such as writing the radial and non-radial parts separately or summing over integer vectors $m$, count, and so does the $xy$ component alone, $P^{\mathcal F}_{xy}=-\frac{\pi^{3/2}}{8V^2}\sum_{k\ne0}\sum_\ell\tilde w_\ell\tilde s_\ell^7e^{-\tilde s_\ell^2|k|^2/4}k_xk_y|\rho(k)|^2$. An expression that omits the $-\tfrac{\tilde s_\ell^2}{2}k\otimes k$ term, uses $\tilde s_\ell^3$ in place of $\tilde s_\ell^5$, or has a different overall prefactor does not count.

### 5

#### Criterion

Identifies the Fourier wavevectors of the triclinic cell, with $h$ the matrix whose columns are the lattice vectors and $m$ integer vectors, as $k=2\pi(h^{-1})^{\top}m$.

#### Category

Scientific Reasoning

#### Weight

2

#### Description

This makes $k\cdot r=2\pi\,m\cdot s$ for fractional coordinates $s$ and reproduces the periodicity of the sheared cell. Equivalent statements count, such as $k=\sum_i m_ib_i$ with reciprocal vectors $a_i\cdot b_j=2\pi\delta_{ij}$, or $k=2\pi h^{-1}m$ when $h$ is defined with the lattice vectors as rows. With the lattice vectors as the columns of $h$, $2\pi h^{-1}m$ gives wrong wavevectors for this non-orthogonal cell.

### 6

#### Criterion

Computes the rescaling factor $\tilde\omega$ of the narrowest retained Gaussian weight as 1.0114 (accepting 1.0110 to 1.0116, or 1.01 when reported to three significant figures).

#### Category

Scientific Reasoning

#### Weight

5

#### Description

Imposing $1/r_c^3=\sum_{\ell=0}^{9}\tilde w_\ell e^{-r_c^2/\tilde s_\ell^2}$ at $r_c=9$ Å, with $\tilde w_0\to\tilde\omega\tilde w_0$, gives 1.011392, which rounds to 1.011 at four and 1.01 at three significant figures. The range also accepts 1.011303, obtained with an untruncated tail, which changes the final answer by less than 0.01 %. A value reported to four or more significant figures outside 1.0110 to 1.0116, such as 1.0106, or $\tilde\omega=1$ (no rescaling), does not count.

### 7

#### Criterion

Computes the short-range (real-space) contribution to $P_{xy}$, summed over all periodic images within $r_c=9$ Å, as $-8.885\times10^{-5}\ e^2/\text{Å}^4$ (within ±0.5 %).

#### Category

Scientific Reasoning

#### Weight

7

#### Description

The rescaled short-range kernel $1/r^3-\sum_\ell\tilde w_\ell e^{-r^2/\tilde s_\ell^2}$ is summed over pair and image separations below $r_c$. Because $r_c$ exceeds half of every perpendicular cell width (9.45, 9.05, 10.20 Å), a minimum-image sum ($\approx-4.79\times10^{-5}$) is wrong. The same value expressed in other pressure units counts.

### 8

#### Criterion

Computes the long-range (Fourier-space) contribution to $P_{xy}$ as $-7.657\times10^{-5}\ e^2/\text{Å}^4$ (within ±0.5 %).

#### Category

Scientific Reasoning

#### Weight

7

#### Description

This is the off-diagonal entry of the converged Fourier-space pressure tensor; integer vectors with $|m_d|\le4$ already converge it, and only the non-radial part contributes. Together with the short-range part it gives $P_{xy}=-1.6542\times10^{-4}\ e^2/\text{Å}^4$. The same value expressed in other pressure units counts. Values outside ±0.5 %, such as $-7.571\times10^{-5}$ from omitting the narrowest-weight rescaling, do not count.

### 9

#### Criterion

States the single-sample variance of an importance-sampled sum $\mu=\sum_kf(k)$, estimated by $f(k)/\mathscr P(k)$ with $k$ drawn from $\mathscr P$, as $\sum_kf(k)^2/\mathscr P(k)-\mu^2$.

#### Category

Scientific Reasoning

#### Weight

5

#### Description

This exact identity, applied with $f$ the $P_{xy}$ summand of the non-radial tensor and $\mathscr P=\mathscr P^{\rm nr}$, gives the requested variance. The identity counts when written for the distribution the response actually samples from, whichever proposal that is (criterion 2 grades the proposal), including equivalent forms such as $\mathbb E_{k\sim\mathscr P}\big[(f(k)/\mathscr P(k))^2\big]-\mu^2$. A denominator that differs from the sampled distribution, the source's variance bound used as an equality, or a missing $-\mu^2$ term does not count.

### 10

#### Criterion

Computes the single-mode variance of the non-radial random-batch estimate of $P_{xy}$ as $5.018\times10^{-7}\ e^4/\text{Å}^8$ (within ±1 %).

#### Category

Scientific Reasoning

#### Weight

8

#### Description

The single-mode estimate is $-S^{\rm nr}|\rho(k)|^2k_xk_y/(8V^2|k|^4)$, whose mean is the long-range $P_{xy}$, where $S^{\rm nr}=\pi^{3/2}\sum_{k\ne0}\sum_\ell\tilde w_\ell\tilde s_\ell^7|k|^4e^{-\tilde s_\ell^2|k|^2/4}$ normalizes $\mathscr P^{\rm nr}$ ($\tilde w_\ell,\tilde s_\ell$ the long-range Gaussian weights, with $\tilde w_0$ rescaled by $\tilde\omega$, and widths). A value of $8.97\times10^{-7}$, from treating the source's variance bound as an equality, or $9.09\times10^{-7}$, from sampling the radial proposal, does not count. The same value expressed in other units of squared pressure counts.

### 11

#### Criterion

Computes the relative standard error of the random-batch estimate of $P_{xy}$ with batch size $P=128$ as 0.3785 (within ±0.5 %, i.e. 0.3766 to 0.3804).

#### Category

Scientific Reasoning

#### Weight

10

#### Description

This equals $\sqrt{5.018\times10^{-7}/128}/|{-1.6542\times10^{-4}}|=6.261\times10^{-5}/1.6542\times10^{-4}$. Omitting the narrowest-weight rescaling gives 0.3748, and treating the variance bound as an equality gives 0.506; neither counts.
