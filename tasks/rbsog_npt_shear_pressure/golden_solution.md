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
