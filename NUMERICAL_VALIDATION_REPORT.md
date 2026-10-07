# Final numerical validation report

This report records the publication-oriented numerical validation for the response-time study. Raw MAT/CSV/PDF artifacts remain attached to the GitHub Actions runs; large raw outputs are not committed to the repository.

## 1. Validated workflow provenance

| Purpose | GitHub Actions run | Commit / branch |
|---|---:|---|
| Refined linear response-time scan | 37466378310 | \`4be6f8eeb98ff0985a686d420982a3f1427bde58\` |
| Isolated-square Hopf/mode selection | 37467535251 | \`b0745936d173612892ccf3d7a7dff65c519a17cb\` |
| Baseline long-time ecology | 37468885257 | \`716f421aee893dd4eddf3fd80361fa6687f03388\` |
| Extended ecology to \(T=1500\) | 37483669030 | \`354865e8567ddaf777ecf75d9988ccabc6d88feb\` |
| Repaired refinement suite | 37483695008 | \`b04496ea18105ad0cd6c26cad985501fae53213b\` |
| Final ecology \(T=3000\) + near-threshold Hopf \(T=4000\) | 37564230418 | \`7ca3d97df31dacec94bccd5863a6a742a812e3c4\` |

The original refinement run 37468890399 and the first final-validation run 37563802808 failed because of workflow command/interpolation errors before producing scientifically meaningful contradictory data. Both were corrected and superseded by successful runs above.

## 2. Linear response-time window

Baseline discretization: cell-centered finite volume, \(100\times100\), \(\Delta t=0.005\), SBDF2/IMEX time stepping. The validation uses tiny pure-mode perturbations and compares measured rates with the *fully discrete* SBDF2 amplification prediction.

| \(\tau\) | measured growth | stability |
|---:|---:|---|
| 0 | -0.185755173 | stable |
| 0.5 | -0.050432152 | stable |
| 0.65 | +0.016731783 | unstable |
| 1 | +0.113451858 | unstable |
| 2 | +0.177617906 | unstable |
| 6 | +0.023763066 | unstable |
| 6.8 | -0.005448571 | stable |
| 8 | -0.043906390 | stable |

Across the eight points, the largest relative discrepancy between measured and fully discrete SBDF2 growth rates is about **0.0123%**, and the largest frequency discrepancy is about **0.000057%**. Hence the numerical calculation cleanly resolves

\[
\text{stable}\;\to\;\text{finite-wavenumber oscillatory instability}\;\to\;\text{restabilization}.
\]

In particular, entry occurs between \(\tau=0.5\) and \(0.65\), while exit occurs between \(\tau=6\) and \(6.8\), in agreement with the analytical continuous-PDE window.

Reliability diagnostics: maximum taxis CFL in the linear validation is about \(1.37\times10^{-2}\), mass-balance residuals are at most \(1.33\times10^{-11}\), and accepted runs remain positive.

## 3. Isolated-square Hopf and mode selection

For a clean nonlinear Hopf test the square is reduced to \(L=5\), preserving the critical physical wavenumber while isolating the first critical family \((2,3)/(3,2)\). The first crossing is

\[
\tau_c \approx 0.608795639,
\]

and the next distinct crossing is delayed to approximately

\[
\tau_{\rm next}\approx0.621178087.
\]

Thus the tests at \(\tau=0.612,0.615,0.618,0.620\) stay inside the isolated first-family regime.

### Cubic normal-form comparison

| branch | \(\tau\) | error in \(|z_1|^2\) | frequency error | field-return residual |
|---|---:|---:|---:|---:|
| pure | 0.612, \(T=4000\) | +18.25% | -0.138% | \(1.70\times10^{-3}\) |
| pure | 0.615, \(T=2000\) | +10.41% | -0.142% | \(1.15\times10^{-3}\) |
| pure | 0.618, \(T=2000\) | +8.34% | -0.147% | \(5.62\times10^{-4}\) |
| pure | 0.620, \(T=2000\) | +7.04% | -0.151% | \(1.21\times10^{-4}\) |
| symmetric | 0.612, \(T=4000\) | +13.29% | -0.141% | \(1.41\times10^{-3}\) |
| symmetric | 0.615, \(T=2000\) | +5.65% | -0.142% | \(5.94\times10^{-4}\) |
| symmetric | 0.618, \(T=2000\) | +3.39% | -0.144% | \(2.76\times10^{-4}\) |
| symmetric | 0.620, \(T=2000\) | +2.01% | -0.145% | \(8.74\times10^{-4}\) |

The nonlinear frequency is therefore predicted very accurately by the cubic normal form (about \(0.15\%\) error), while the amplitude is **semi-quantitative rather than a precision match**. The amplitude errors are consistent with finite-distance/higher-order effects and slow near-threshold relaxation; they should not be described as sub-percent agreement.

### Transverse behavior at \(\tau=0.612\)

For the pure branch the cubic prediction gives a negative transverse rate (\(-3.07\times10^{-5}\)). The measured transverse contamination remains at only a few \(10^{-6}\) and the late log-slope is negative, consistent with transverse decay.

For the symmetric branch the cubic prediction gives a positive transverse rate (\(+2.81\times10^{-4}\)). In the \(T=4000\) run, the transverse ratio increases to about \(2.0\times10^{-5}\); its late log-slope is positive (order \(5\times10^{-4}\)). Thus the **sign and order of the transverse tendency are resolved**, but the measured exponent should not be advertised as a high-precision normal-form coefficient.

Hopf reliability diagnostics remain clean: maximum CFL is below \(7.3\times10^{-3}\), mass-balance residuals are below \(8\times10^{-13}\), and accepted runs remain positive.

## 4. Long-time ecological effects

The publication baseline uses \(L=10\), \(100\times100\), \(\Delta t=0.005\), upwind taxis flux. The final \(\tau=1,2\) calculations were extended to \(T=3000\). The reported values below are late-time averages on \(t\in[2250,3000]\).

The homogeneous reference (\(\tau=0\), and again after restabilization at \(\tau=8\)) is approximately

\[
C_0=24.14903843,\qquad
P_{\rm tot,0}=53.97259494,\qquad
N_{\rm tot,0}=40.77524218,\qquad
O_0=1.
\]

| \(\tau\) | consumption | \(P_{\rm tot}\) | \(N_{\rm tot}\) | overlap \(O\) |
|---:|---:|---:|---:|---:|
| 0 | 24.149038 | 53.972595 | 40.775242 | 1.000000 |
| 1 | 23.565178 | 52.543833 | 41.053638 | 0.971182 |
| 2 | 23.342497 | 51.762068 | 41.192292 | 0.956336 |
| 8 | 24.149038 | 53.972595 | 40.775242 | 1.000000 |

Relative to the homogeneous reference:

- \(\tau=1\): mean consumption changes by **-2.418%**, predator total by **-2.647%**, and prey total by **+0.683%**.
- \(\tau=2\): mean consumption changes by **-3.340%**, predator total by **-4.096%**, and prey total by **+1.023%**.

### Exact predation decomposition

For
\[
C=|\Omega|\bar P F(\bar N)+|\Omega|\bar P\{\overline{F(N)}-F(\bar N)\}
+|\Omega|{\rm Cov}_x(P,F(N)),
\]
the decomposition closes to roundoff.

At \(\tau=1\),
\[
\Delta C=-0.583861
=(-0.522123)_{\rm abundance}
+(-0.174771)_{\rm heterogeneity}
+(+0.113033)_{\rm association}.
\]

At \(\tau=2\),
\[
\Delta C=-0.806542
=(-0.816306)_{\rm abundance}
+(-0.238430)_{\rm heterogeneity}
+(+0.248194)_{\rm association}.
\]

Hence the reduction in mean predation is **not** caused by a negative spatial association term. In these runs, abundance change and Holling-II heterogeneity/saturation reduce consumption, while the positive spatial-association contribution partially offsets that reduction.

### Late-time stability

Comparing \(t\in[2000,2500]\) with \(t\in[2500,3000]\):

- \(\tau=1\): consumption, \(P_{\rm tot}\), \(N_{\rm tot}\), and overlap differ by at most about **0.11%**.
- \(\tau=2\): the same primary quantities differ by at most about **0.19%**.

Pattern amplitudes retain modest residual modulation (roughly \(1.6\%-2.5\%\) between those windows), and cue mismatch changes by about \(1.9\%\) at \(\tau=1\) and \(3.7\%\) at \(\tau=2\). Therefore primary integrated ecological quantities are numerically settled to well below \(1\%\), while pattern amplitude and cue-mismatch values should be reported as oscillatory diagnostics rather than sub-percent constants.

## 5. Numerical refinement

### Threshold sign robustness

For \(\tau=0.65\) the baseline measured growth is \(+0.01673170\). With \(\Delta t/2\) it is \(+0.01660779\), and with \(150^2\) cells it is \(+0.01674434\). The sign is unchanged.

For \(\tau=6.8\) the baseline measured growth is \(-0.00544859\). With \(\Delta t/2\) it is \(-0.00547806\), and with \(150^2\) cells it is \(-0.00543879\). The sign is unchanged.

The magnitude change is below about \(0.75\%\) for these growth rates; frequency changes remain below about \(0.13\%\).

### Nonlinear ecological refinement at \(\tau=1\)

At \(T=300\):

- halving \(\Delta t\) changes consumption by about **0.019%**, \(N_{\rm tot}\) by **0.002%**, and overlap by **0.022%**;
- increasing from \(100^2\) to \(150^2\) changes consumption by about **0.468%**, \(N_{\rm tot}\) by **0.132%**, and overlap by **0.569%**.

The core ecological effects therefore satisfy the requested \(\sim1\%\) publication-level consistency criterion.

Cue mismatch changes by about **8.3%** under spatial refinement. It is retained only as a secondary diagnostic and should not carry a precision claim.

## 6. Numerical reliability

Across the accepted publication-oriented runs:

- no silent negative-value clipping is used;
- accepted states remain nonnegative;
- taxis CFL values remain comfortably below the diagnostic safety threshold;
- discrete mass balances close at roughly \(10^{-11}\) or better in the large-domain runs and about \(10^{-13}\) in the isolated Hopf tests;
- solver/workflow failures were treated as computational failures, never as PDE blow-up evidence.

No continuation/hysteresis calculation is required for the present manuscript claims: the completed calculations do not provide affirmative evidence that multiple attractors materially affect the stated response-time-window, local Hopf, or integrated ecological conclusions.

## 7. Manuscript-safe numerical statements

The following statements are supported by the completed calculations:

1. The numerical method resolves the predicted finite response-time window, including both entry and restabilization.
2. Linear modal growth rates and oscillation frequencies agree essentially exactly with the fully discrete SBDF2 prediction at the tested points.
3. On the isolated \(L=5\) square, the first Hopf family produces small-amplitude periodic branches whose nonlinear frequencies agree with the cubic normal form to about \(0.15\%\).
4. Cubic amplitude predictions capture the scale and branch trend but are only semi-quantitative (roughly \(2\%-18\%\) errors over the tested points).
5. Pure-branch transverse perturbations decay, while the symmetric branch exhibits the predicted positive transverse tendency; this is a qualitative stability validation, not a high-precision exponent measurement.
6. In the unstable response-time regime, long-time integrated ecological means show reduced mean consumption and predator abundance, increased prey abundance, and reduced normalized predator-prey overlap.
7. The exact predation decomposition shows that abundance and Holling-II heterogeneity/saturation effects drive the reduction in consumption, while positive spatial association partially compensates it.
8. The main ecological quantities pass time-step and spatial-grid refinement at the \(\lesssim1\%\) level; cue mismatch is more mesh sensitive and should be treated as secondary.

## 8. Recommended final figures

- **Figure 1:** analytical response-time growth envelope + eight measured linear validation points, with entry/exit thresholds.
- **Figure 2:** isolated-square pure/symmetric Hopf branches: measured vs cubic amplitude scale, frequency comparison, and one transverse-stability diagnostic.
- **Figure 3:** representative \(P,N,W\) states/dynamics for \(\tau=0,1,2,8\).
- **Figure 4:** long-time consumption / prey total / overlap together with the three-term predation decomposition.
- **Supplement:** refinement tables, CFL/mass/positivity diagnostics, long histories, and cue-mismatch mesh sensitivity.

## 9. Bottom line

The numerical section is sufficiently validated for the present SIAP manuscript claims. The strongest quantitative claims are the response-time window validation, nonlinear frequency validation, long-time integrated ecological effects, exact predation decomposition, and refinement of core observables. The main caveats are that cubic branch amplitudes are semi-quantitative rather than precision matches and cue mismatch is spatial-grid sensitive.
