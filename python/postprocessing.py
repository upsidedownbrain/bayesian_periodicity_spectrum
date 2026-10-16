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

"""Posterior-odds transforms matching the MATLAB implementation."""
import numpy as np


def posterior_odds_smoothing(bf_psi, f, KernelSigma=0, EpsClip=1e-12):
    bf_psi = np.maximum(np.asarray(bf_psi, dtype=float), np.log10(EpsClip))
    f = np.asarray(f, dtype=float).ravel()
    if KernelSigma > 0:
        D = np.abs(f[:, None]-f[None, :])
        K = np.exp(-.5*(D/KernelSigma)**2)
        K /= np.maximum(np.sum(K, axis=1, keepdims=True), np.finfo(float).eps)
        return K@bf_psi
    return bf_psi


def posterior_odds_temperature(bf_psi, EpsClip=1e-12, Temp=1):
    bf_psi = np.maximum(np.asarray(bf_psi, dtype=float), np.log10(EpsClip))
    vector = bf_psi.ndim == 1
    if vector: bf_psi = bf_psi[:, None]
    F, M = bf_psi.shape
    post = np.empty((F, M))
    total = np.sum(bf_psi, axis=1)
    for i in range(M):
        prior = (total-bf_psi[:, i])/max(M-1, 1)
        post[:, i] = bf_psi[:, i]+Temp*prior
    return post[:, 0] if vector else post
