# ============================
# bayesian_periodicity package
# ============================

# ----------------------------
# core.py
# ----------------------------
import numpy as np
from scipy.signal import stft
from scipy.signal import correlate
import pywt

from .statistics import ttest_bf_directional, credible_interval_individual_spectrum


def periodicity_analysis(s, Fs, window_BF, method, psi_non_osc,
                         window=None, overlapLength=None, FFTLength=None,
                         VoicesPerOctave=16, CI=False, BF_corr=True):

    s = np.asarray(s).flatten()

    if window is None:
        window = Fs
    if overlapLength is None:
        overlapLength = window - 1
    if FFTLength is None:
        FFTLength = Fs

    if method == 'STFT':
        f, t, Zxx = stft(s, fs=Fs, window='hann', nperseg=window,
                         noverlap=overlapLength, nfft=FFTLength, boundary=None)
        pos = (f > 0) & (f < Fs/2)
        freqs = f[pos]
        tfr = Zxx[pos, :]

    elif method == 'CWT':
        widths = np.geomspace(1, 128, VoicesPerOctave * 6)
        tfr, freqs = pywt.cwt(s, widths, 'morl', sampling_period=1/Fs)
        mask = freqs >= 1
        tfr = tfr[mask, :]
        freqs = freqs[mask]
    else:
        raise ValueError('Method must be STFT or CWT')

    phi = np.angle(tfr)
    dphi = np.diff(phi, axis=1)

    omega0 = 2 * np.pi * freqs[:, None] / Fs
    resid = dphi - omega0
    dphi_corr = (resid + np.pi) % (2*np.pi) - np.pi

    omega_inst = omega0 + dphi_corr
    f_inst = (Fs/(2*np.pi)) * omega_inst
    f_inst = f_inst.T

    segments = f_inst.shape[0] // window_BF
    power = []

    for z in range(segments):
        seg = tfr[:, z*window_BF:(z+1)*window_BF]
        power.append(np.mean(np.abs(seg)**2, axis=1))

    power = np.array(power).T

    psi = []
    for z in range(segments):
        seg = f_inst[z*window_BF:(z+1)*window_BF, :]
        psi.append(1 / np.std(seg, axis=0))

    psi = np.array(psi).T
    psi_norm = psi / psi_non_osc['mean_psi_non_osc_fitted'][:, None]

    results = {}
    results['freqs'] = freqs
    results['power_sp'] = np.mean(power, axis=1)
    results['psi_norm'] = np.mean(psi_norm, axis=1)

    if CI:
        mu_mean, mu_ci, df = credible_interval_individual_spectrum(psi_norm.T, freqs)
        results['mu_mean'] = mu_mean
        results['mu_ci'] = mu_ci
        results['df'] = df

    psi_shift = psi_norm - 1
    bf_psi = []

    for z in range(len(freqs)):
        data = psi_shift[z, :]

        if BF_corr:
            ac = correlate(data - np.mean(data), data - np.mean(data), mode='full')
            rho1 = ac[len(data)] / np.max(ac)
    return results