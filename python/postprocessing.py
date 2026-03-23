# ============================
# bayesian_periodicity package
# ============================

# ----------------------------
# postprocessing.py
# ----------------------------
import numpy as np


def posterior_odds_smoothing(bf_psi, f, KernelSigma=0, EpsClip=1e-12):
    bf_psi = np.maximum(bf_psi, np.log10(EpsClip))

    if KernelSigma > 0:
        D = np.abs(f[:, None] - f[None, :])
        K = np.exp(-0.5 * (D / KernelSigma)**2)
        K = K / np.maximum(np.sum(K, axis=1, keepdims=True), np.finfo(float).eps)
        return K @ bf_psi
    else:
        return bf_psi


def posterior_odds_temperature(bf_psi, EpsClip=1e-12, Temp=1):
    bf_psi = np.maximum(bf_psi, np.log10(EpsClip))

    F, M = bf_psi.shape
    postOdds = np.zeros_like(bf_psi)

    sumLogs = np.sum(bf_psi, axis=1)

    for i in range(M):
        logPrior = (sumLogs - bf_psi[:, i]) / max(M - 1, 1)
        logPrior = Temp * logPrior
        postOdds[:, i] = logPrior + bf_psi[:, i]

    return postOdds