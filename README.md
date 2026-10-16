# Bayesian periodicity spectrum

This repository contains MATLAB and Python implementations of the **Bayesian periodicity spectrum**, a method to quantify oscillatory activity based on the variability of the instantaneous frequency. The Bayesian periodicity spectrum provides a frequency-resolved measure of oscillatory activity that is robust to confounds from non-oscillatory components. Additionally, it incorporates a Bayesian framework to formally assess the **presence or absence of oscillations** at each frequency. This approach addresses key limitations of traditional methods that rely solely on the power spectrum for separating oscillatory and non-oscillatory components.

---

# Method overview

The Bayesian periodicity spectrum is computed through the following steps:

1. Time-frequency decomposition (STFT or CWT)  
2. Instantaneous frequency estimation (phase-based)  
3. Variability-based periodicity metric  
4. Normalization using non-oscillatory activity  
5. Bayesian inference using Bayes factors  

---

# Repository structure
matlab/

├── periodicity_analysis.m

├── non_oscillatory_periodicity_spectrum.m

├── credible_interval_individual_spectrum.m

├── ttest_bf.m

├── spectra_comparison.m

├── posterior_odds_smoothing.m

├── posterior_odds_temperature.m

python/

├── core.py

├── statistics.py

├── postprocessing.py

├── simulation.py

├── __init__.py

├── pyproject.toml

└── tests/

---

# Installation

Clone this repository, then install the Python package from its root:

```shell
python -m pip install ./python
```

For MATLAB, add the `matlab` folder to the MATLAB path. See
[the Python documentation](python/README.md) for API and array conventions.

# Basic usage
Given a signal of interest s(t), the "Bayesian periodicity spectrum toolbox" first calculates the corresponding periodicity of non-oscilaltory components, depending of the method used for the frequency component decomposition — for this toolbox, the short-time Fourier transform (STFT) or the continuous wavelet transform (CWT). If multiple signals with different durations share the same decomposition parameters, use the shortest signal duration to compute the periodicity of the non-oscillatory components. This reduces computation time — shorter simulated signals require fewer resources — and the resulting periodicity can then be applied to all signals of interest. The periodicity spectrum is estimated as the inverse square root of the variance of the instantaneous frequency for each component. Instantaneous frequency traces are divided into segments of equal length in order to estimate the Bayesian periodicity spectrum. This strategy allows to assess the level of evidence for the presence or the absence of oscillatory activity in a single signal.

Below there is a typical workflow to obtain the power, periodicity, and Bayesian periodicity spectra of a signal in MATLAB. This workflow is equivalent in Python.

## Example (MATLAB)

```matlab

% Signal
Fs = 128;
window_BF = 128;
method = 'CWT';

% --- Example signal ---
t = 0:1/Fs:5-1/Fs;
s = sin(2*pi*10*t) + 0.5*randn(size(t));

% --- Step 1: non-oscillatory model ---
psi_non_osc = non_oscillatory_periodicity_spectrum(length(s), Fs, window_BF, method);

% --- Step 2: periodicity analysis ---
results = periodicity_analysis(s, Fs, window_BF, method, psi_non_osc);

% Outputs
power_sp = results.power_sp;
psi_norm = results.psi_norm;
bf_psi   = results.bf_psi;
freqs    = results.freqs;

```

## Example (Python)

```python
import numpy as np
from bayesian_periodicity import periodicity_analysis, non_oscillatory_periodicity_spectrum

# Signal
Fs = 128
window_BF = 128
method = "CWT"

# Example signal
t = np.arange(0, 5, 1/Fs)
s = np.sin(2*np.pi*10*t) + 0.5*np.random.randn(len(t))

# --- Step 1: non-oscillatory model ---
psi_non_osc = non_oscillatory_periodicity_spectrum(len(s), Fs, window_BF, method)

# --- Step 2: periodicity analysis ---
results = periodicity_analysis(s, Fs, window_BF, method, psi_non_osc)

power_sp = results["power_sp"]
psi_norm = results["psi_norm"]
bf_psi   = results["bf_psi"]
freqs    = results["freqs"]

```

To compare periodicity (or power/Bayesian spectra) between two populations:
- MATLAB: spectra_comparison.m
- Python: spectra_comparison (in statistics.py)
This function performs a Bayesian paired t-test at each frequency component.

---

# Notes on MATLAB vs Python
The MATLAB implementation serves as the reference implementation. Python follows its STFT conventions, analytic Morlet CWT, segment sample standard deviations, directional JZS Bayes factors, ESS correction, and AR(1) credible intervals. The supplied regression tests compare results against MATLAB R2021a outputs within floating-point tolerances. Independently generated Monte Carlo references can differ because the languages use different random-number streams.

`bf_psi` contains log10 Bayes factors; `bf10_psi` contains raw Bayes factors. `window_BF` is measured in original signal samples, including when the STFT hop is greater than one. Use matching decomposition options when simulating a reference and analysing a signal. Credible intervals are conditional on the estimated AR(1) correlation; constant columns are rejected because their posterior under the reference prior is improper.

MATLAB requirements:
- Signal Processing Toolbox
- Wavelet Toolbox
- Statistics and Machine Learning Toolbox
- DSP System Toolbox
- Curve Fitting Toolbox (for the CWT simulation reference fit)

Python requirements:
- numpy>=1.23
- scipy>=1.9

PyWavelets is no longer required. Python must be version 3.8 or newer.

## Verification

From the repository root:

```shell
python python/tests/test_equivalence.py
```

The eight regression groups cover STFT/CWT pipelines, nonunit hops, all tails, extreme t statistics, ESS, missing data, AR(1) intervals, paired comparisons, posterior odds, and simulation fits. The MATLAB reference outputs are bundled; running the Python tests does not require MATLAB. To regenerate the fixtures with MATLAB and the toolboxes listed above, run `python/tests/export_matlab_reference.m`.

---

# License
MIT License

---

# Reference
Please, if you use this code, cite this publication:

---
