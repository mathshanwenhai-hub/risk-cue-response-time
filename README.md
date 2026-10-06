# Risk-cue response time in a predator-prey taxis system

Reproducible MATLAB code for a study of **response-time-induced oscillatory instability** in a predator-prey pursuit-evasion model with a predator-generated risk cue.

## Model

The numerical experiments use the three-field PDE system

- predator diffusion and prey tracking,
- prey diffusion and avoidance of the predator-generated cue,
- a dynamic cue with response-time parameter \(\tau\),

with homogeneous Neumann boundary conditions and the Holling-II response
\(F(N)=aN/(1+ahN)\).

The benchmark domain is \([0,10]^2\).

## Benchmark prediction

For the parameter set in `matlab/riskcue_parameters.m`, the continuous PDE mode calculation predicts approximately

- lower response-time boundary: `tau_in = 0.6087956390`;
- upper response-time boundary: `tau_out = 6.6445349603`;
- first critical square-domain mode family: `(4,6)` and `(6,4)`;
- last exiting mode family: `(4,8)` and `(8,4)`.

These values are re-computed by the CI theory workflow rather than hard-coded into the simulation.

## Numerical method

The nonlinear solver is a cell-centered finite-volume discretization with

- implicit diffusion;
- explicit taxis and reaction terms;
- SBDF2 / IMEX-BDF2 after an Euler startup step;
- no silent clipping of negative values;
- diagnostics for taxis CFL, discrete mass balances, minimum values, modal amplitudes, cue mismatch, overlap, and predation decomposition.

The baseline publication configuration is \(N_x=N_y=100\) and \(\Delta t=0.005\).

The upwind taxis flux is used for finite-amplitude simulations. Centered flux is used only in targeted near-bifurcation tests where numerical upwind diffusion would contaminate the cubic amplitude comparison.

## Reproduce locally

With MATLAB on the path:

```matlab
addpath('matlab');
github_theory_assertions;
github_smoke_test;
```

A single production parameter point can be run with

```matlab
github_run_case(1,'window',100,0.005,180,'upwind',0);
```

Outputs are written to `results/` and `figures/`.

## GitHub Actions

The repository contains the following workflows:

1. **01 Theory and smoke test** — runs automatically after numerical-code changes.
2. **02 Response-time window scan** — manual 8-point scan over `tau = 0, 0.5, 0.65, 1, 2, 6, 6.8, 8`.
3. **03 Hopf branch and mode selection** — manual pure/symmetric near-entry runs with a transverse seed.
4. **04 Ecological diagnostics** — manual long-time runs at `tau = 0, 1, 2, 8`.
5. **05 Numerical refinement** — manual time-step and spatial-grid checks at `tau = 0.65, 1, 6.8`.
6. **06 Optional continuation** — manual sequential response-time continuation.

Each production job uploads MAT files, CSV histories/summaries, logs, and generated PDF figures as GitHub Actions artifacts.

## Important interpretation notes

- A stopped numerical run is **not** interpreted as PDE blow-up.
- Long-time averages are provisional until adjacent time windows and refinement runs agree.
- Linear growth-rate and frequency fits are diagnostics and are compared against the discrete FV spectrum.
- The repository contains numerical/reproducibility code only; the manuscript and private research notes are intentionally not included.

## MATLAB version

GitHub Actions are pinned to MATLAB **R2025b** for reproducibility. The source code itself avoids specialized MATLAB toolboxes.

## License

No open-source license has been assigned yet. All rights are reserved until a license is added explicitly.
