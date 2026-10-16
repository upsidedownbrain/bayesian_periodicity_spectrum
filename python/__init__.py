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

"""Bayesian periodicity spectrum, MATLAB-compatible Python implementation."""
from .core import periodicity_analysis
from .simulation import non_oscillatory_periodicity_spectrum
from .statistics import (ttest_bf, ttest_bf_directional,
                         credible_interval_individual_spectrum,
                         spectra_comparison, PrecisionNotAchieved)
from .postprocessing import posterior_odds_smoothing, posterior_odds_temperature

__version__ = '2.0.0'
__all__ = ['periodicity_analysis', 'non_oscillatory_periodicity_spectrum',
           'ttest_bf', 'ttest_bf_directional', 'credible_interval_individual_spectrum',
           'spectra_comparison', 'posterior_odds_smoothing',
           'posterior_odds_temperature', 'PrecisionNotAchieved']
