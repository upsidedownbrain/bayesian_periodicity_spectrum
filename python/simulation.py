# Authors: Jesús Pardo-Valencia, Guglielmo Foffani
# Version: 2.0.0 
# Date: October 2026
# Repository: GitHub
# License: MIT License
# Copyright (c) 2025 jpv
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Non-oscillatory references using the same transforms as the analysis."""
import numpy as np

from .core import _decompose, _instantaneous_frequency, _segments, _periodicity, _integer
from .statistics import _real


def non_oscillatory_periodicity_spectrum(signal_duration, Fs, window_BF, method,
                                         window=None, hop=1, FFTLength=None,
                                         VoicesPerOctave=16, N_simulations=1000,
                                         *, random_state=None):
    """Monte Carlo reference and STFT constant/CWT log-poly6 fit.

    random_state accepts a NumPy generator or seed for reproducible Python
    runs. MATLAB and NumPy random generators are not expected to agree.
    A CWT fit is returned as descending polynomial coefficients (fit_mean).
    """
    durations = _real(signal_duration, 'signal_duration').ravel()
    if durations.size not in (1, 2) or np.any(durations <= 0):
        raise ValueError('signal_duration requires one positive length or [min,max].')
    if durations.size == 2 and durations[0] > durations[1]:
        raise ValueError('The first duration must not exceed the second.')
    lengths = np.floor(np.linspace(durations[0], durations[-1], 10)+.5).astype(int)
    nsim = _integer(N_simulations, 'N_simulations')
    if nsim % 10:
        raise ValueError('N_simulations must be a multiple of ten, as required by the MATLAB reference.')
    hop = _integer(hop, 'hop')
    window = int(np.floor(Fs+.5)) if window is None else window
    FFTLength = int(np.floor(Fs+.5)) if FFTLength is None else FFTLength
    rng = np.random.default_rng(random_state)
    models, common_freqs = [], None
    for length in lengths:
        accum = []
        for _ in range(nsim//10):
            s = rng.standard_normal(length)
            s = (s-np.mean(s))/np.std(s, ddof=1)
            tfr, freqs, actual_hop = _decompose(s, Fs, method, window, window-hop,
                                                FFTLength, VoicesPerOctave)
            inst = _instantaneous_frequency(tfr, freqs, Fs, actual_hop)
            edges = _segments(inst, window_BF, actual_hop, minimum=1)
            accum.append(np.mean(_periodicity(inst, edges), axis=1))
        if common_freqs is None:
            common_freqs = freqs
        if len(freqs) < len(common_freqs) or not np.allclose(freqs[:len(common_freqs)], common_freqs, rtol=1e-12, atol=1e-12):
            raise ValueError('Simulation lengths do not share the reference frequency grid.')
        models.append(np.mean(accum, axis=0)[:len(common_freqs)])
    mean = np.mean(models, axis=0)
    if not np.all(np.isfinite(mean)) or np.any(mean <= 0):
        raise ValueError('Simulation produced an invalid periodicity reference.')
    if method == 'STFT':
        if len(mean) < 3: raise ValueError('STFT reference fit requires at least three positive frequency bins.')
        fit = np.mean(mean[1:-1])
        fitted = np.full(mean.shape, fit)
    else:
        if len(mean) < 7: raise ValueError('CWT poly6 fit requires at least seven frequencies.')
        fit = np.polyfit(np.log(common_freqs), np.log(mean), 6)
        fitted = np.exp(np.polyval(fit, np.log(common_freqs)))
    return dict(freqs=common_freqs, mean_psi_non_osc=mean, fit_mean=fit,
                mean_psi_non_osc_fitted=fitted)
