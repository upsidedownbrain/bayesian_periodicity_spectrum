# Bayesian periodicity spectrum — Python 2.0.0

This implementation follows the files in the repository's
[`matlab` directory](../matlab) (October 2026).
It uses NumPy and SciPy; PyWavelets is no longer required. The analytic
Morlet transform uses MATLAB's `amor` frequency grid, L1 normalization and
reflection extension, rather than the real `morl` wavelet.

## Installation

```powershell
python -m pip install ./python
```

```python
import numpy as np
from bayesian_periodicity import (
    periodicity_analysis, non_oscillatory_periodicity_spectrum, ttest_bf,
    credible_interval_individual_spectrum,
)

rng = np.random.default_rng(12)
s = rng.normal(size=6000)
reference = non_oscillatory_periodicity_spectrum(
    len(s), 100, 100, "STFT", N_simulations=1000, random_state=12,
)
result = periodicity_analysis(s, 100, 100, "STFT", reference, CI=True)
bf, p, t, n, logbf, diagnostics = ttest_bf(
    T=2, N=20, rho1=.2, tail="right", return_diagnostics=True,
)
```

## API and numerical conventions

- `ttest_bf` returns `(bf10, pValue, T, N)` by default. `return_log=True`
  adds `log10bf`; `return_diagnostics=True` adds both `log10bf` and the
  diagnostic dictionary. MATLAB-style name/value calls such as
  `ttest_bf('T', 2, 'N', 20, 'tail', 'right')` also work.
- `ttest_bf_directional` remains available with its four original outputs
  and default `tail='right'`. `ttest_bf` defaults to `'both'`, like MATLAB.
- One-sample, paired and pooled-variance independent tests omit NaNs.
  ESS is supported for one-sample/paired tests and updates T, N and df
  together, with the same bounds as MATLAB. Independent direct T uses
  `N=[n1,n2]` and scalar T.
- `Dim=1/2` and `SampleDimension=1/2` use MATLAB numbering, not zero-based
  NumPy axes. Vectors of outputs are 1-D NumPy arrays. `psi_segments` and
  `power_segments` have shape `(frequency, segment)`; `f_inst` has shape
  `(time, frequency)`; `mu_ci` has shape `(frequency, 2)`.
- Credible intervals return `(mu_mean, mu_ci, df)` for compatibility;
  `return_rho=True` adds rho1. They use the AR(1) GLS Student-t posterior
  conditional on estimated rho1. Constant columns are rejected because
  the reference-prior posterior is improper. Square inputs need an
  explicit `SampleDimension`.
- STFT uses an unscaled FFT and a periodic Hann window, without boundary
  padding or an incomplete final frame. Instantaneous frequency accounts
  for the hop. `window_BF` is in original signal samples, and segment
  boundaries are rounded up to the next available frame. Each segment
  requires at least two frames. Use identical transform options for the
  simulated reference and the analysed signal.
- CWT power is divided by frequency. Segment standard deviations use
  `ddof=1`. All MATLAB result fields, including `df` and `df_ci` when CI is
  enabled, are provided. `bf_psi` stores log10 BF, not raw BF.
- Simulation defaults to 1000 trials and ten length groups. As with the
  MATLAB reference, `N_simulations` must be divisible by ten. The STFT
  reference is a constant fitted from the interior frequency bins; CWT
  uses degree-six regression of log periodicity on log frequency.
  `fit_mean` is a scalar for STFT and seven descending polynomial
  coefficients for CWT (instead of MATLAB's `cfit` object).
- NumPy and MATLAB use different random-number streams: independently
  generated Monte Carlo reference values will differ. `random_state`
  permits reproducible Python simulations. Identical supplied noise was
  used in the cross-language tests to verify the simulation calculations.
- Posterior-odds routines retain MATLAB's current clipping and averaging
  rules, including their behavior for nonfinite inputs.

Cross-language verification uses MATLAB R2021a reference outputs for
STFT/CWT pipelines, odd FFT and signal lengths, nonunit/nondividing hops,
all tails, extreme t statistics, missing data, ESS, AR(1) intervals,
paired comparisons, simulation fits and posterior-odds calculations.
Floating-point equality is tested with tolerances, not bitwise identity.

Run the supplied eight regression groups from the source folder:

```powershell
python tests/test_equivalence.py
```

Run the installation command from the repository root and the test command
from the `python` source folder.

The MATLAB outputs used by those tests are included in
`tests/python_reference.mat`, so running the Python tests does not require
MATLAB. `tests/export_matlab_reference.m` records how the fixture was
generated. Run it with MATLAB to regenerate the reference file in `tests`.
