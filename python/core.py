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

"""Periodicity analysis with MATLAB-compatible transform conventions."""
import numpy as np
from scipy.signal.windows import hann

from .statistics import (_real, _lag1, credible_interval_individual_spectrum,
                         ttest_bf)


def _integer(value, name, minimum=1):
    if not np.isscalar(value) or not np.isfinite(value) or value != int(value) or value < minimum:
        raise ValueError(name + ' must be an integer >= ' + str(minimum))
    return int(value)


def _decompose(s, Fs, method, window, overlapLength, FFTLength, VoicesPerOctave):
    if not np.isscalar(Fs) or not np.isfinite(Fs) or Fs <= 0:
        raise ValueError('Fs must be finite and positive.')
    if method not in ('STFT', 'CWT'): raise ValueError('method must be STFT or CWT.')
    s = _real(s, 's').ravel(order='F')
    if method == 'STFT':
        window = _integer(Fs if window is None else window, 'window')
        overlap = _integer(Fs-1 if overlapLength is None else overlapLength, 'overlapLength', 0)
        nfft = _integer(Fs if FFTLength is None else FFTLength, 'FFTLength')
        hop = window-overlap
        if hop <= 0 or window > len(s) or nfft < window:
            raise ValueError('Require overlapLength < window <= signal length and FFTLength >= window.')
        frames = np.lib.stride_tricks.sliding_window_view(s, window)[::hop]
        spectrum = np.fft.fft(frames*hann(window, sym=False), n=nfft, axis=1).T
        freqs = np.fft.fftfreq(nfft, 1/Fs)
        keep = (freqs > 0) & (freqs < Fs/2)
        return spectrum[keep], freqs[keep], hop
    voices = _integer(VoicesPerOctave, 'VoicesPerOctave')
    if not 1 <= voices <= 48: raise ValueError('VoicesPerOctave must be between 1 and 48.')
    if len(s) < 4: raise ValueError('CWT requires at least four signal samples.')
    # Analytic Morlet: positive-frequency Gaussian centred at omega=6,
    # L1 normalization (peak 2), 10% Nyquist cutoff, reflection extension.
    smallest = (6+np.sqrt(2*np.log(10)))/np.pi
    largest = len(s)/(2*np.sqrt(2))
    count = int(np.floor(voices*np.log2(largest/smallest)))+1
    scales = smallest*2**(np.arange(count)/voices)
    freqs = 6*Fs/(2*np.pi*scales)
    keep = freqs >= 1
    scales, freqs = scales[keep], freqs[keep]
    pad = len(s)//2
    extended = np.pad(s, (pad, pad), mode='symmetric')
    omega = 2*np.pi*np.fft.fftfreq(len(extended))
    if len(extended) % 2 == 0: omega[len(extended)//2] = np.pi
    transformed = np.fft.fft(extended)
    tfr = np.empty((len(freqs), len(s)), dtype=complex)
    chunk = max(1, 1000000//len(extended))
    for start in range(0, len(freqs), chunk):
        a = scales[start:start+chunk, None]*omega
        response = 2*np.exp(-.5*(a-6)**2)*(omega > 0)
        tfr[start:start+chunk] = np.fft.ifft(response*transformed, axis=1)[:, pad:pad+len(s)]
    return tfr, freqs, 1


def _instantaneous_frequency(tfr, freqs, Fs, hop):
    phase_difference = np.diff(np.angle(tfr), axis=1)
    omega = 2*np.pi*freqs[:, None]/Fs
    residual = (phase_difference-omega*hop+np.pi) % (2*np.pi)-np.pi
    return ((omega+residual/hop)*Fs/(2*np.pi)).T


def _segments(f_inst, window_BF, hop, minimum=3):
    width = _integer(window_BF, 'window_BF')
    count = int(np.floor(len(f_inst)*hop/width))
    if count < minimum: raise ValueError('At least %d complete segments are required.' % minimum)
    edges = np.ceil(np.arange(count+1)*width/hop).astype(int)
    if np.any(np.diff(edges) < 2):
        raise ValueError('Each segment requires at least two instantaneous-frequency frames. Increase window_BF or decrease the STFT hop.')
    return edges


def _periodicity(f_inst, edges):
    with np.errstate(divide='ignore'):
        return np.stack([1/np.std(f_inst[a:b], axis=0, ddof=1)
                         for a, b in zip(edges[:-1], edges[1:])], axis=1)


def periodicity_analysis(s, Fs, window_BF, method, psi_non_osc,
                         window=None, overlapLength=None, FFTLength=None,
                         VoicesPerOctave=16, CI=False, BF_corr=True):
    """MATLAB-equivalent spectra; vectors are 1-D NumPy arrays.

    psi_segments/power_segments are frequency x segment; f_inst is time x
    frequency; mu_ci is frequency x 2. bf_psi contains log10 directional BF.
    """
    tfr, freqs, hop = _decompose(s, Fs, method, window, overlapLength, FFTLength, VoicesPerOctave)
    f_inst = _instantaneous_frequency(tfr, freqs, Fs, hop)
    edges = _segments(f_inst, window_BF, hop)
    base = _real(psi_non_osc['mean_psi_non_osc_fitted'], 'mean_psi_non_osc_fitted').ravel()
    if base.size != len(freqs) or np.any(base <= 0):
        raise ValueError('Reference must contain one positive finite value per analysed frequency.')
    power = np.stack([np.mean(np.abs(tfr[:, a:b])**2, axis=1)
                      for a, b in zip(edges[:-1], edges[1:])], axis=1)
    if method == 'CWT': power = power/freqs[:, None]
    psi = _periodicity(f_inst, edges)/base[:, None]
    result = dict(freqs=freqs, tfr=tfr, f_inst=f_inst, n_segments=len(edges)-1,
                  n_timepoints_used=int(edges[-1]), power_segments=power,
                  power_sp=np.mean(power, axis=1), psi_segments=psi,
                  psi_norm=np.mean(psi, axis=1))
    if CI:
        mu, ci, df, rho = credible_interval_individual_spectrum(psi, freqs, SampleDimension=2, return_rho=True)
        result.update(mu_mean=mu, mu_ci=ci, df=df, df_ci=df, rho1_ci=rho)
    rho = np.array([_lag1(row-1) if BF_corr else 0.0 for row in psi])
    bf, p, t, n, logbf = ttest_bf(psi-1, Dim=2, rho1=rho if BF_corr else None,
                                 tail='right', return_log=True)
    result.update(bf_psi=logbf, bf10_psi=bf, rho1=rho, T_bf=t,
                  N_eff_bf=n, p_bf=p)
    return result
