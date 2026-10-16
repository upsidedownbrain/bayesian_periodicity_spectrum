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

"""Bayesian tests and AR(1) intervals matching the MATLAB implementation."""
from functools import lru_cache
import warnings

import numpy as np
from scipy import integrate, optimize, special, stats


class PrecisionNotAchieved(RuntimeError):
    """The Bayes-factor integral did not meet the requested tolerance."""


def _real(value, name, finite=True):
    a = np.asarray(value)
    if a.dtype.kind not in 'iuf' or np.iscomplexobj(a) or not a.size:
        raise ValueError(name + ' must contain real numeric values.')
    a = a.astype(float)
    if finite and not np.all(np.isfinite(a)):
        raise ValueError(name + ' must be finite.')
    return a


def _shape_parameter(a, shape, name):
    a = _real(a, name)
    if a.size == 1:
        return a.reshape(())
    if a.size != np.prod(shape):
        raise ValueError(name + ' must be scalar or have one value per test.')
    return a.reshape(shape, order='F')


def ttest_bf(X=None, *args, **kwargs):
    """Return (BF10, p, T, N), with MATLAB sample dimensions (Dim=1/2).

    Use ttest_bf(T=..., N=...) or ttest_bf('T', ..., 'N', ...).
    return_log=True adds log10BF; return_diagnostics=True adds both log10BF
    and the diagnostic dictionary. The first four outputs retain their order.
    rho1 applies ESS to nominal one-sample/paired T and N, exactly once.
    """
    cfg = dict(tail='both', scale=np.sqrt(.5), T=None, N=None, rho1=None,
               paired=True, nQuad=16, RelTol=1e-10, MaxQuad=None, Dim=None,
               return_log=False, return_diagnostics=False)
    aliases = {k.lower(): k for k in cfg}
    Y = 0
    if isinstance(X, str):
        args = (X,) + args
        X = None
    elif args and not isinstance(args[0], str):
        Y, args = args[0], args[1:]
    if len(args) % 2:
        raise ValueError('Options must be name/value pairs.')
    supplied = list(zip(args[::2], args[1::2])) + list(kwargs.items())
    seen = set()
    for key, value in supplied:
        canonical = aliases.get(str(key).lower())
        if canonical is None or canonical in seen:
            raise ValueError('Unknown or duplicate option: ' + str(key))
        cfg[canonical] = value
        seen.add(canonical)
    tail = str(cfg['tail']).lower()
    if tail not in ('right', 'left', 'both'):
        raise ValueError('tail must be right, left, or both.')
    r, tol = float(cfg['scale']), float(cfg['RelTol'])
    if not np.isfinite(r) or r <= 0 or not 1e-12 <= tol <= 1e-3:
        raise ValueError('Invalid scale or RelTol.')
    q0 = cfg['nQuad']
    qmax = cfg['MaxQuad'] if cfg['MaxQuad'] is not None else max(128, 2*q0)
    if not (8 <= q0 <= 128 and int(q0) == q0 and 16 <= qmax <= 256
            and int(qmax) == qmax and qmax > q0):
        raise ValueError('Invalid nQuad/MaxQuad.')
    if not isinstance(cfg['paired'], (bool, np.bool_)):
        raise ValueError('paired must be boolean.')
    dim = cfg['Dim']
    if dim is not None and dim not in (1, 2):
        raise ValueError('Dim must be 1 or 2.')
    if cfg['T'] is not None:
        if X is not None or dim is not None:
            raise ValueError('Direct T cannot be combined with data or Dim.')
        T = _real(cfg['T'], 'T')
        N = _real(cfg['N'], 'N')
        if T.ndim > 2 or np.any(N <= 0):
            raise ValueError('Invalid T or N.')
        if N.size == 2 and T.size == 1:
            if np.any(N < 2) or np.any(N != np.floor(N)):
                raise ValueError('Independent sample sizes must be integers >= 2.')
            n1, n2 = N.ravel()
            df, N, use_ess = n1+n2-2, n1*n2/(n1+n2), False
        else:
            if np.any(N <= 1):
                raise ValueError('N must be greater than one.')
            df, use_ess = N-1, True
    else:
        if cfg['N'] is not None:
            raise ValueError('N is calculated from data; use it only with direct T.')
        X, Y = _real(X, 'X', False), _real(Y, 'Y', False)
        if X.ndim > 2 or Y.ndim > 2 or np.any(np.isinf(X)) or np.any(np.isinf(Y)):
            raise ValueError('Data must be at most 2-D and cannot contain Inf.')
        if dim is None and (X.ndim < 2 or 1 in X.shape) and (Y.size == 1 or Y.ndim < 2 or 1 in Y.shape):
            X = X.reshape(-1, 1, order='F')
            if Y.size != 1:
                Y = Y.reshape(-1, 1, order='F')
        else:
            if X.ndim < 2: X = np.atleast_2d(X)
            if Y.size != 1 and Y.ndim < 2: Y = np.atleast_2d(Y)
        axis = dim-1 if dim is not None else next((i for i, n in enumerate(X.shape) if n != 1), 0)
        if Y.size == 1:
            if not np.isfinite(Y.item()): raise ValueError('Null mean must be finite.')
            values = X-Y.item()
            use_ess = True
        elif cfg['paired']:
            if X.shape != Y.shape: raise ValueError('Paired arrays must have identical shapes.')
            values, use_ess = X-Y, True
        else:
            if X.shape[1-axis] != Y.shape[1-axis]: raise ValueError('Independent arrays must match outside Dim.')
            n1, n2 = np.sum(~np.isnan(X), axis=axis), np.sum(~np.isnan(Y), axis=axis)
            if np.any(n1 < 2) or np.any(n2 < 2): raise ValueError('At least two observations per sample are required.')
            df = n1+n2-2
            N = n1*n2/(n1+n2)
            pooled = ((n1-1)*np.nanvar(X, axis=axis, ddof=1)+(n2-1)*np.nanvar(Y, axis=axis, ddof=1))/df
            T = (np.nanmean(X, axis=axis)-np.nanmean(Y, axis=axis))/np.sqrt(pooled/N)
            use_ess = False
        if use_ess:
            N = np.sum(~np.isnan(values), axis=axis)
            if np.any(N < 2): raise ValueError('At least two valid observations/pairs are required.')
            df = N-1
            with np.errstate(divide='ignore', invalid='ignore'):
                T = np.nanmean(values, axis=axis)/np.sqrt(np.nanvar(values, axis=axis, ddof=1)/N)
        T = np.asarray(T)
        # Python vectors remain 1-D; a single vector test returns scalars.
        if T.size == 1: T = T.reshape(())
    T = np.asarray(T, dtype=float)
    if not np.all(np.isfinite(T)): raise ValueError('DegenerateData: t statistic is not finite.')
    shape = T.shape
    N = _shape_parameter(N, shape, 'N')
    df = _shape_parameter(df, shape, 'df')
    if cfg['rho1'] is not None:
        if not use_ess: raise ValueError('rho1 is supported only for one-sample/paired tests.')
        rho = _shape_parameter(cfg['rho1'], shape, 'rho1')
        if np.any(np.abs(rho) >= 1) or np.any(N < 2): raise ValueError('Invalid rho1 or nominal N.')
        rho = np.clip(rho, -.99, .99)
        neff = np.maximum(2, np.minimum(N, N*(1-rho)/(1+rho)))
        T = T*np.sqrt(neff/N)
        N, df = neff, neff-1
    if tail == 'both': p = 2*stats.t.sf(np.abs(T), df)
    elif tail == 'right': p = stats.t.sf(T, df)
    else: p = stats.t.cdf(T, df)
    tv = T.ravel()
    nv, dv = np.broadcast_to(N, shape).ravel(), np.broadcast_to(df, shape).ravel()
    logs, diagnostics = _backend(tv, nv, dv, r, tail, int(q0), tol, int(qmax))
    with np.errstate(over='ignore', under='ignore'):
        bf = np.exp(logs).reshape(shape)
    out = (bf, np.asarray(p).reshape(shape), T, N)
    if cfg['return_log'] or cfg['return_diagnostics']:
        out += (logs.reshape(shape)/np.log(10),)
    if cfg['return_diagnostics']:
        out += ({key: val.reshape(shape) for key, val in diagnostics.items()},)
    return out


def ttest_bf_directional(X=None, tail='right', scale=np.sqrt(.5), rho1=None, **kwargs):
    """Backward compatible four-output directional test."""
    return ttest_bf(X, tail=tail, scale=scale, rho1=rho1, **kwargs)


def _log_beta_half(a):
    a = np.asarray(a)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        approx = .5*np.log(np.pi)-.5*np.log(a)+1/(8*a)-1/(192*a**3)+1/(640*a**5)
    return np.where(a > 100, approx, special.betaln(a, .5))


def _clamp(a):
    return np.where(np.abs(a) < 1e-300, np.where(a < 0, -1e-300, 1e-300), a)


def _log_beta_cf(a, x):
    qab, qap, qam = a+.5, a+1, a-1
    c = np.ones_like(a)
    d = 1/_clamp(1-qab*x/qap)
    h = d.copy()
    active = np.ones(a.shape, dtype=bool)
    for m in range(1, 4097):
        ix = np.flatnonzero(active)
        if not ix.size: break
        ai, xi, m2 = a[ix], x[ix], 2*m
        aa = m*(.5-m)*xi/((qam[ix]+m2)*(ai+m2))
        di, ci = 1/_clamp(1+aa*d[ix]), _clamp(1+aa/c[ix])
        hi = h[ix]*di*ci
        aa = -(ai+m)*(qab[ix]+m)*xi/((ai+m2)*(qap[ix]+m2))
        di, ci = 1/_clamp(1+aa*di), _clamp(1+aa/ci)
        delta = di*ci
        h[ix], d[ix], c[ix] = hi*delta, di, ci
        active[ix[np.abs(delta-1) < 5e-14]] = False
    if np.any(active) or np.any(h <= 0) or not np.all(np.isfinite(h)):
        raise PrecisionNotAchieved('Directional incomplete beta did not converge.')
    return np.log(h)


def _log_integrand(v, t, n, nu, r, tail):
    v = np.atleast_1d(v)
    t, n, nu = (np.atleast_1d(a).reshape(-1, 1) for a in (t, n, nu))
    loga = np.log(n)+2*np.log(r)+v
    log1pa = np.logaddexp(0, loga)
    with np.errstate(divide='ignore'):
        logt = np.log(np.abs(t))
    logratio = 2*logt-np.log(nu)
    logd0, logdg = np.logaddexp(0, logratio), np.logaddexp(0, logratio-log1pa)
    logx = logdg-logd0
    alpha = (nu+1)/2
    base = -.5*log1pa-alpha*logx
    if tail != 'both':
        if tail == 'left': t = -t
        logz = logt+.5*(loga-log1pa)+.5*(np.log1p(1/nu)-logdg)
        z = np.exp(np.minimum(logz, np.log(np.finfo(float).max)))
        prob = stats.t.cdf(z*np.sign(t), nu+1)
        with np.errstate(divide='ignore'):
            ydir = base+np.log(2)+np.log(prob)
        mask = (t < 0) & ((-alpha*logx > 32) | (prob == 0))
        if np.any(mask):
            aa = np.broadcast_to(alpha, base.shape)[mask]
            xx = np.exp(logx[mask])
            cf = np.zeros_like(xx)
            nz = xx > 0
            cf[nz] = _log_beta_cf(aa[nz], xx[nz])
            ydir[mask] = (-.5*log1pa[mask]+.5*np.log1p(-xx)
                          -_log_beta_half(aa)-np.log(aa)+cf)
        base = ydir
    return base-.5*np.log(2*np.pi)-.5*v-.5*np.exp(-v)


def _tail_bounds(t, n, nu, r, tail, lo, hi):
    with np.errstate(divide='ignore'):
        c = (nu+1)/2*np.logaddexp(0, 2*np.log(np.abs(t))-np.log(nu))
    factor = 0
    if tail != 'both':
        factor = np.log(2)
        c = np.where((t if tail == 'right' else -t) <= 0, 0, c)
    upper = c-.5*(np.log(n)+2*np.log(r))-.5*np.log(2*np.pi)-hi+factor
    lower = c+.5*lo-.5*np.exp(-lo)+.5*np.log(2/np.pi)+factor
    return lower, upper


@lru_cache(maxsize=16)
def _legendre(q):
    return np.polynomial.legendre.leggauss(q)


def _quadrature(t, n, nu, r, tail, q, lo, hi):
    x, w = _legendre(q)
    edges = np.arange(lo, hi, 2)
    v = (edges[:, None]+1+x).ravel()
    lw = np.tile(np.log(w), len(edges))
    result = np.empty(len(t))
    chunk = max(1, 1000000//len(v))
    for start in range(0, len(t), chunk):
        ix = slice(start, start+chunk)
        result[ix] = special.logsumexp(_log_integrand(v, t[ix], n[ix], nu[ix], r, tail)+lw, axis=1)
    return result


def _adaptive(t, n, nu, r, tail, tol):
    lo = -12
    hi = max(48, 2*np.ceil((2*np.log(max(abs(t), np.finfo(float).tiny))-np.log(n)-2*np.log(r)+40)/2))
    def fun(v): return _log_integrand(v, t, n, nu, r, tail).ravel()
    for _ in range(12):
        scan = np.linspace(lo, hi, max(61, int(np.ceil((hi-lo)/2))+1))
        i = int(np.argmax(fun(scan)))
        opt = optimize.minimize_scalar(lambda v: -fun(v)[0], bounds=(scan[max(0, i-1)], scan[min(len(scan)-1, i+1)]), method='bounded')
        mode, peak = opt.x, -opt.fun
        points = np.unique(np.r_[mode+np.arange(-8, 9, 2), (lo+hi)/2])
        points = points[(points > lo) & (points < hi)]
        f = lambda v: np.exp(fun(v)[0]-peak)
        v1, e1 = integrate.quad(f, lo, hi, points=points, epsabs=tol/32, epsrel=max(tol/32, 100*np.finfo(float).eps), limit=1000)
        v2, e2 = integrate.quad(f, lo, hi, points=points, epsabs=tol/128, epsrel=max(tol/128, 100*np.finfo(float).eps), limit=1000)
        if v2 > 0:
            value = peak+np.log(v2)
            lb, ub = _tail_bounds(t, n, nu, r, tail, lo, hi)
            err = abs(v1/v2-1)+max(e1, e2)/v2+np.exp(lb-value)+np.exp(ub-value)
            if np.isfinite(value) and err <= tol: return value, err
        lo, hi = lo-2, hi+16
    raise PrecisionNotAchieved('Adaptive JZS integration did not converge.')


def _backend(t, n, nu, r, tail, q0, tol, qmax):
    size = len(t)
    logs, errors = np.full(size, np.nan), np.full(size, np.inf)
    nodes, fallback = np.zeros(size, dtype=int), np.zeros(size, dtype=bool)
    lo = -12
    hi = max(48, 2*np.ceil(np.max(2*np.log(np.maximum(np.abs(t), np.finfo(float).tiny))-np.log(n)-2*np.log(r)+40)/2))
    pending = np.arange(size)
    if hi-lo <= 2048:
        for _ in range(8):
            prev = _quadrature(t[pending], n[pending], nu[pending], r, tail, q0, lo, hi)
            q = q0
            while q < qmax and pending.size:
                q = min(2*q, qmax)
                cur = _quadrature(t[pending], n[pending], nu[pending], r, tail, q, lo, hi)
                with np.errstate(over='ignore'):
                    e = np.abs(np.expm1(prev-cur))
                lb, ub = _tail_bounds(t[pending], n[pending], nu[pending], r, tail, lo, hi)
                trunc = np.exp(lb-cur)+np.exp(ub-cur)
                ok = np.isfinite(cur) & (e <= tol/4) & (trunc <= tol/4)
                logs[pending[ok]], errors[pending[ok]], nodes[pending[ok]] = cur[ok], e[ok]+trunc[ok], q
                pending, prev = pending[~ok], cur[~ok]
            if not pending.size: break
            lb, ub = _tail_bounds(t[pending], n[pending], nu[pending], r, tail, lo, hi)
            if np.all(np.exp(lb-prev)+np.exp(ub-prev) <= tol/4): break
            lo, hi = lo-2, hi+16
            if hi-lo > 2048: break
    for j in pending:
        logs[j], errors[j] = _adaptive(t[j], n[j], nu[j], r, tail, tol)
        fallback[j] = True
    if not np.all(np.isfinite(logs)) or np.any(errors > tol):
        raise PrecisionNotAchieved('BF integration did not meet RelTol.')
    return logs, dict(converged=np.ones(size, dtype=bool), relativeErrorEstimate=errors,
                      nodesPerPanel=nodes, adaptiveFallback=fallback,
                      bfOverflow=logs > np.log(np.finfo(float).max),
                      bfUnderflow=logs < np.log(np.finfo(float).tiny))


def _lag1(x):
    centered = x-np.mean(x)
    if np.std(centered, ddof=1) == 0: return 0.0
    rho = np.dot(centered[:-1], centered[1:])/np.dot(centered, centered)
    return float(rho) if np.isfinite(rho) else 0.0


def credible_interval_individual_spectrum(X, freqs, alpha=.05, *, SampleDimension=None,
                                          Alpha=None, return_rho=False):
    """AR(1) GLS Student-t interval conditional on estimated rho1.

    SampleDimension uses MATLAB numbering: 1=segments x frequencies,
    2=frequencies x segments. Square inputs require an explicit dimension.
    Returns the original (mean, CI, df) tuple; return_rho adds rho1.
    """
    if Alpha is not None: alpha = Alpha
    if not 0 < alpha < 1: raise ValueError('Alpha must be between zero and one.')
    X, freqs = _real(X, 'X'), _real(freqs, 'freqs').ravel()
    if X.ndim == 1: X = X[:, None]
    if X.ndim != 2: raise ValueError('X must be 2-D.')
    f = len(freqs)
    if SampleDimension is None:
        matches = [X.shape[1] == f, X.shape[0] == f]
        if all(matches): raise ValueError('AmbiguousOrientation: specify SampleDimension for square inputs.')
        if not any(matches): raise ValueError('Frequency vector length does not match X.')
        SampleDimension = 1 if matches[0] else 2
    if SampleDimension not in (1, 2) or X.shape[2-SampleDimension] != f:
        raise ValueError('Frequency size does not match SampleDimension.')
    if SampleDimension == 2: X = X.T
    n, f = X.shape
    if n < 3: raise ValueError('At least 3 segments are required.')
    means, ci, rho = np.empty(f), np.empty((f, 2)), np.empty(f)
    df = n-1
    one = np.ones(n)
    distance = np.abs(np.arange(n)[:, None]-np.arange(n))
    for j, x in enumerate(X.T):
        if np.all(x == x[0]): raise ValueError('DegenerateData: constant frequency column; reference-prior posterior is improper.')
        rho[j] = np.clip(_lag1(x), -.99, .99)
        R = rho[j]**distance
        information = one@np.linalg.solve(R, one)
        mu = (one@np.linalg.solve(R, x))/information
        residual = x-mu
        sse = residual@np.linalg.solve(R, residual)
        if not np.isfinite(sse) or sse <= 0: raise ValueError('DegenerateData: residual sum of squares must be positive and finite.')
        means[j] = mu
        width = stats.t.ppf(1-alpha/2, df)*np.sqrt(sse/df/information)
        ci[j] = [mu-width, mu+width]
    result = (means, ci, df)
    return result+(rho,) if return_rho else result


def spectra_comparison(spectra_1, spectra_2, tail='both'):
    """Paired Bayesian comparisons, one row per frequency."""
    x, y = _real(spectra_1, 'spectra_1', False), _real(spectra_2, 'spectra_2', False)
    if x.ndim != 2 or x.shape != y.shape: raise ValueError('Spectra must be equally sized frequency x sample matrices.')
    return np.array([ttest_bf(a, b, tail=tail)[0] for a, b in zip(x, y)])
