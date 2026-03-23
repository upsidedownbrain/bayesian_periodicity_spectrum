# ============================
# bayesian_periodicity package
# ============================

# ----------------------------
# simulation.py
# ----------------------------
import numpy as np
from .core import periodicity_analysis


def non_oscillatory_periodicity_spectrum(signal_duration, Fs, window_BF, method,
                                          window=None, hop=1, FFTLength=None,
                                          VoicesPerOctave=16, N_simulations=10000):

    if isinstance(signal_duration, (int, float)):
        lengths = [int(signal_duration)] * 10
    elif len(signal_duration) == 2:
        lengths = np.linspace(signal_duration[0], signal_duration[1], 10).astype(int)
    else:
        raise ValueError('signal_duration must have 1 or 2 elements')

    psi_all = []

    for L in lengths:
        sims = np.random.randn(L, N_simulations // len(lengths))
        sims = (sims - sims.mean(axis=0)) / sims.std(axis=0)

        for i in range(sims.shape[1]):
            res = periodicity_analysis(sims[:, i], Fs, window_BF, method,
                                       {'mean_psi_non_osc_fitted': np.ones(100)})
            psi_all.append(res['psi_norm'])

    psi_all = np.array(psi_all)

    mean_psi = np.mean(psi_all, axis=0)

    results = {}
    results['mean_psi_non_osc'] = mean_psi
    results['mean_psi_non_osc_fitted'] = mean_psi

    return results
