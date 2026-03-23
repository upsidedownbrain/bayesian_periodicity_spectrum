# ============================
# bayesian_periodicity package
# ============================

# ----------------------------
# statistics.py
# ----------------------------
import numpy as np
            df = Neff - 1
            N = Neff
    else:
        df = N - 1

    if tail == 'both':
        p_value = 2 * (1 - t.cdf(abs(T), df))
    elif tail == 'right':
        p_value = 1 - t.cdf(T, df)
    else:
        p_value = t.cdf(T, df)

    def t_likelihood(delta):
        return (1 + (T - np.sqrt(N) * delta)**2 / df) ** (-(df + 1) / 2)

    def half_cauchy(delta):
        return (2 / np.pi / scale) * (1 / (1 + (delta / scale)**2))

    integrand = lambda d: t_likelihood(d) * half_cauchy(d)

    ml_H1, _ = quad(integrand, 0, np.inf, epsabs=1e-10, epsrel=1e-6)
    ml_H0 = t_likelihood(0)

    bf10 = ml_H1 / ml_H0

    return bf10, p_value, T, N


def credible_interval_individual_spectrum(X, freqs, alpha=0.05):
    X = np.asarray(X)
    freqs = np.asarray(freqs)

    if X.shape[1] != len(freqs):
        X = X.T

    E = X.shape[0]
    if E < 3:
        raise ValueError("Need at least 3 epochs.")

    mu_mean = np.mean(X, axis=0)
    s = np.std(X, axis=0, ddof=1)
    se = s / np.sqrt(E)

    df = E - 1
    tcrit = t.ppf(1 - alpha / 2, df)

    ci_low = mu_mean - tcrit * se
    ci_high = mu_mean + tcrit * se

    mu_ci = np.vstack((ci_low, ci_high)).T

    return mu_mean, mu_ci, df


def spectra_comparison(spectra_1, spectra_2, tail='both'):
    spectra_1 = np.asarray(spectra_1)
    spectra_2 = np.asarray(spectra_2)

    bf10 = []
    for j in range(spectra_1.shape[0]):
        bf, _, _, _ = ttest_bf_directional(spectra_1[j, :] - spectra_2[j, :], tail=tail)
        bf10.append(bf)

    return np.array(bf10)