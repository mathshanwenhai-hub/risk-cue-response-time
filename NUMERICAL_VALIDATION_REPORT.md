# Numerical validation report

## Validation status

Publication-oriented computations have been checked from GitHub Actions artifacts.

Validated runs:
- linear window: 37466378310 (4be6f8eeb98ff0985a686d420982a3f1427bde58)
- isolated L=5 Hopf branches: 37467535251 (b0745936d173612892ccf3d7a7dff65c519a17cb)
- baseline ecology: 37468885257 (716f421aee893dd4eddf3fd80361fa6687f03388)
- extended ecology: 37483669030 (354865e8567ddaf777ecf75d9988ccabc6d88feb)
- repaired refinement: 37483695008 (b04496ea18105ad0cd6c26cad985501fae53213b)

Run 37468890399 is superseded: it failed before meaningful computation because of a workflow command-formatting error.

## Linear response-time window

The measured stability sequence is stable at tau=0 and 0.5, unstable at tau=0.65, 1, 2 and 6, and stable again at tau=6.8 and 8. Across the eight validation points, the largest relative discrepancy between measured and fully discrete SBDF2 growth rates is 0.0123%, and the largest frequency discrepancy is 0.00006%. Maximum taxis CFL is 0.0137, maximum mass-balance residual is 1.33e-11, and all recorded minima remain positive.

## Isolated-square Hopf validation

On L=5 the first critical family (2,3)/(3,2), tau_c approximately 0.60879564, is separated from the next crossing near 0.62117809. For tau=0.612, 0.615, 0.618 and 0.620, nonlinear-frequency errors relative to the cubic normal form are about 0.135-0.151%. Amplitude-squared errors are about 13-15% closest to onset and decrease to about 2-7% farther from threshold. Full-field return residuals are 1.21e-4 to 1.15e-3.

Pure-branch transverse contamination decreases. Symmetric-branch transverse growth is qualitatively visible for tau at least 0.615; tau=0.612 is too close to threshold for a precise measured transverse exponent. Maximum Hopf CFL is 0.00725 and maximum mass residual is 7.94e-13.

## Long-time ecological means

| tau | consumption | Ptotal | Ntotal | overlap |
|---:|---:|---:|---:|---:|
| 0 | 24.149038 | 53.972595 | 40.775242 | 1.000000 |
| 1 | 23.568088 | 52.545333 | 41.056089 | 0.971102 |
| 2 | 23.334715 | 51.743083 | 41.195859 | 0.955955 |
| 8 | 24.149038 | 53.972595 | 40.775242 | 1.000000 |

Relative to tau=0, tau=1 changes consumption by -2.406%, Ptotal by -2.644%, and Ntotal by +0.689%. At tau=2 the corresponding changes are -3.372%, -4.131%, and +1.032%.

The exact predation decomposition closes to roundoff. At tau=1, total consumption change -0.580951 equals abundance -0.520395 plus heterogeneity -0.174075 plus spatial association +0.113519. At tau=2, -0.814323 equals -0.823407 -0.240625 +0.249708. Thus spatial association partly offsets the consumption decrease rather than causing it.

Between late windows 1000-1250 and 1250-1500, consumption, Ptotal, Ntotal and overlap differ by less than 0.4% for tau=1 and 2. Pattern amplitudes and cue mismatch are more window-sensitive: the largest changes are about 1.5% at tau=1 and 5.8% at tau=2. These secondary diagnostics should not be reported as sub-percent converged constants on the present horizon.

## Refinement

For tau=0.65 and 6.8, both dt/2 and N=150 checks preserve the growth-rate sign. Growth-rate changes are below 0.75% and frequency changes below 0.13%.

For nonlinear tau=1, dt/2 changes consumption, Ntotal and overlap by 0.019%, 0.002% and 0.022%. N=150 changes them by 0.468%, 0.132% and 0.569%. Cue mismatch changes by about 8.3% under spatial refinement and should be treated as a secondary mesh-sensitive diagnostic.

## Manuscript-safe statements

1. The computation resolves a response-time instability window with entry between tau=0.5 and 0.65 and restabilization between tau=6 and 6.8.
2. Measured linear growth rates and frequencies agree extremely closely with the fully discrete SBDF2 prediction.
3. The isolated L=5 Hopf family has nonlinear frequency within about 0.15% of cubic normal-form prediction and field-return residuals of order 1e-4 to 1e-3.
4. Cubic amplitude scaling is reproduced semi-quantitatively and improves away from onset.
5. Pure transverse perturbations decay; symmetric transverse growth is a qualitative observation, not a measured exponent.
6. At tau=1 and 2, integrated late ecological means show lower consumption and predator abundance but higher prey abundance relative to the homogeneous reference.
7. Abundance and heterogeneity reduce consumption while positive spatial association partially offsets the reduction.
8. Core ecological observables pass the requested refinement criterion; cue mismatch is more mesh sensitive.
9. Accepted runs remain positive, have small CFL values, and preserve discrete mass balances to about 1e-11 or better.

## Figure plan

Figure 1: response-time growth envelope and eight linear validation points.
Figure 2: L=5 pure/symmetric Hopf amplitude and frequency comparisons.
Figure 3: representative tau=0,1,2,8 dynamics or snapshots.
Figure 4: late ecological means plus predation decomposition.
Supplement: refinement, mass/CFL/positivity diagnostics, long histories, and cue-mismatch sensitivity.

## Remaining caveat

The primary ecological observables meet the late-window criterion at T=1500. A stricter requirement that every non-negligible pattern diagnostic change by less than about 1% would require a longer tau=1/2 ecology run. The current execution environment did not permit launching that longer workflow. This does not alter the integrated ecological conclusions above.
