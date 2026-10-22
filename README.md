# fractional-semilinear-calderon

Matlab code for the paper "A Numerical Method for the Fractional Semilinear Calderón Problem" by Mukul Dwivedi, Jesse Railo, and Andreas Rupp.

## Running the code

Unzip to a new directory, make it current in MATLAB R2025b, and enter:

```matlab
run_all
```

No toolbox is required. This archive contains MATLAB source and this README only. It includes no measurements, saved solutions, figures or replay script.

Both examples for $p=1,2,3$ are computed and saved before any figure export starts. Figures remain visible and open.

## Numerical problem

$$(-\Delta)^s u + q u^p = 0$$

in the unit interval/disk, with exterior input $f$. $q$ is a reaction coefficient.

- $d=1,2$; $s=0.35$; $p=1,2,3$.
- Known even symmetry in 1D and radial symmetry in 2D are imposed throughout. The 2D operator and physical weights are two-dimensional; angular variation of an unknown coefficient is outside these benchmarks.
- $b(t)=\exp\bigl(1-1/(1-t^2)\bigr)$ for $|t|<1$, and zero otherwise.
- $q(r)=b(r/0.75)$; $f(r)=16\,b\bigl((r-1.6)/0.2\bigr)$; reconstruction radius $D=0.8$.
- Observation radii: 1.05 to 1.3 (120 cells), 2 to 3 (80 cells).
- $h$ is radial shell width, not a general 2D finite-element diameter.
- $h=0.04, 0.02, 0.01, 0.005, 0.0025$; $\delta=10^{-1}, 10^{-3}, 10^{-5}, 10^{-7}, 10^{-9}, 10^{-10}$.
- At $h=0.0025$: 400 radial state cells, 320 coefficient cells, 200 observations; $B$ is 200-by-400. The forward basis has Jacobi indices 0 through 256: 257 coefficients, polynomial degree 256 in $r^2$ before the boundary weight.

## Noise and parameter rules

Noise has physical weighted norm $\delta\,\|g\|$, where $g$ is the full projected noiseless response. The implementation uses $\mu$ + noise, algebraically equivalent to $(g+\text{noise})-g_0$; it does not add artificial subtraction roundoff. $\delta$ is relative; the absolute data error also includes forward/preprocessing error. Five fixed Gaussian directions are paired across noise levels and meshes. No seed or regularization parameter is selected using true errors.

- **Noise study:** $h=0.0025$, $\tau=2.5\cdot10^{-7}$, largest $\alpha$ meeting weighted discrepancy $1.05\,\delta\,\|g\|$, $\beta=10^{-10}\,\mathrm{mean}_D\bigl(|\hat U|^{2p}\bigr)$.
- **Mesh study:** $\delta=10^{-10}$; each direction's finest-grid $\alpha$, $\tau$ and $\beta$ are held fixed across all meshes. Coarse-grid discrepancies can exceed threshold.

All 30 noise and 25 mesh runs per dimension and power are retained: 330 runs. The same five directions are also paired across $p=1,2,3$. A failed noise criterion stops loudly; completed inverse records are saved in `partial_results.mat` with `complete=false` if a later inverse call fails.

Exact $q$, $u$ and reaction arrays generate data and evaluate errors; after data generation they do not affect reconstruction or parameter selection. The known-input lift is not diagnostic truth. $q_H$ is unconstrained and may be negative. The state multiplier is clipped to $[-16,16]$; weak differentiation uses the unclipped state.

## Interpretation

This is a smooth radial benchmark with known geometry/order/symmetry, radial noise, accurately known input and a prescribed synthetic noise norm. It does not demonstrate general angular recovery or the theorem's simultaneous limit. The observation space stays fixed at 200 cells and $\tau$ is fixed in both studies. No logarithmic exponent is fitted or claimed for the algorithm. A coefficient fit on $D$ is not a bound for the nonlinear PDE residual on all of $\Omega$. Native MATLAB execution/export verification remains necessary.
