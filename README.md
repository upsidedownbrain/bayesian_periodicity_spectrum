# Bayesian Periodicity Spectrum

This repository contains MATLAB and Python implementations of the **Bayesian periodicity spectrum**, a method to quantify oscillatory activity based on the variability of the instantaneous frequency. The Bayesian periodicity spectrum provides a frequency-resolved measure of oscillatory activity that is robust to confounds from non-oscillatory components. Additionally, it incorporates a Bayesian framework to formally assess the **presence or absence of oscillations** at each frequency. This approach addresses key limitations of traditional methods that rely solely on the power spectrum for separating oscillatory and non-oscillatory components.

---

# Method Overview

The Bayesian periodicity spectrum is computed through the following steps:

1. Time-frequency decomposition (STFT or CWT)  
2. Instantaneous frequency estimation (phase-based)  
3. Variability-based periodicity metric  
4. Normalization using non-oscillatory activity  
5. Bayesian inference using Bayes factors  

---

# Repository Structure
matlab/

├── periodicity_analysis.m
├── non_oscillatory_periodicity_spectrum.m
├── credible_interval_individual_spectrum.m
├── ttest_bf_directional.m
├── spectra_comparison.m
├── posterior_odds_smoothing.m
├── posterior_odds_temperature.m

python/
├── core.py
├── statistics.py
├── postprocessing.py
├── simulation.py
├── init.py

---

# Basic usage
Given a signal of interest s(t), the "Bayesian periodicity spectrum toolbox" first calculates the corresponding periodicity of non-oscilaltory components, depending of the method used for the frequency component decomposition — for this toolbox, the short-time Fourier transform (STFT) or the continuous wavelet transform (CWT). If multiple signals with different durations share the same decomposition parameters, use the shortest signal duration to compute the periodicity of the non-oscillatory components. This reduces computation time — shorter simulated signals require fewer resources — and the resulting periodicity can then be applied to all signals of interest. The periodicity spectrum is estimated as the inverse square root of the variance of the instantaneous frequency for each component. Instantaneous frequency traces are divided into segments of equal length in order to estimate the Bayesian periodicity spectrum. This strategy allows to assess the level of evidence for the presence or the absence of oscillatory activity in a single signal.

Below there is a typical workflow to obtain the power, periodicity, and Bayesian periodicity spectra of a signal in MATLAB. This workflow is equivalent in Python.

## Example (MATLAB)

```matlab
% ---------------------------------------------------------------------
% In this example the signal is stored in a variable called 's', 
% a M x 1 array.
% The sampling frequency is set to 128 Hz, the segments in which the
% time series will be divided is set to a length of 128 samples, and
% the decomposition method is the continuous wavelet transform.
% Optional input arguments are set to default values.
% ---------------------------------------------------------------------
% Signal
Fs = 128;
window_BF = 128;
method = 'CWT';

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
The MATLAB implementation serves as the reference implementation. The Python version is scientifically equivalent, but minor numerical differences may occur due to differences from STFT/CWT implementations.

---

# License
MIT License

---

# Reference
Please, if you use this code, cite this publication:
Pardo-Valencia J., Ammann C., Fernández-García C., Alonso-Frech F., Foffani G. (2026).  
*The Bayesian periodicity spectrum to quantify oscillatory activity in biological signals.*  
Science Advances.

---
