function [postOdds] = posterior_odds_smoothing(bf_psi, f, varargin)
% This function estimates posterior odds from Bayes factors using a spatial
% smoothing prior over frequency bins.
%
% Authors: Jesús Pardo-Valencia, Guglielmo Foffani
% Version: 2.0.0
% Date: March 2026
% Reference: Pardo-Valencia et al. (2026). The Bayesian periodicity
% spectrum to quantify oscillatory activity in biological signals. Science
% Advances.
% Repository: GitHub
% License: MIT License
% Copyright (c) 2025 jpv
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
% bf_psi: [F x M] (frecuencies x signals) in log-scale
% f: [F x 1] in Hz
% Optional:
p = inputParser;
p.addParameter('KernelSigma', 0);   % in Hz; 0 = no smoothing
p.addParameter('EpsClip', 1e-12);   % avoids -Inf/NaN
p.parse(varargin{:});
sigma = p.Results.KernelSigma;
epsc  = p.Results.EpsClip;

[F,M] = size(bf_psi);
bf_psi = max(bf_psi, log10(epsc));   % avoids -Inf
postOdds = zeros(F,M);

% --- Smoothing (optional) ---
if sigma>0
    % gaussian kernel on frequencies
    D = abs(f(:) - f(:)');
    K = exp(-0.5*(D/sigma).^2);
    K = K ./ max(sum(K,2), eps);
    postOdds = K * bf_psi;
else
    postOdds = bf_psi;
end
end