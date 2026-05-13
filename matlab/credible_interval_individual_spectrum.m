function [mu_mean, mu_ci, df] = credible_interval_individual_spectrum(X, freqs, varargin)
% bayes_ci_individual_spectrum
%
% Bayesian credible intervals for the periodicity spectrum
% at the INDIVIDUAL level (across epochs).
%
% INPUT
%   X     : [epochs x frequencies] periodicity values
%   freqs : frequency vector
%
% OPTIONS
%   'Alpha' (0.05)
%
% OUTPUT
%   .freqs
%   .mu_mean   [1 x F] posterior mean
%   .mu_ci     [F x 2] credible intervals
%   .df        degrees of freedom (epochs-1)
%
% Authors: Jesús Pardo-Valencia, Guglielmo Foffani
% Version: 2.0.0
% Date: March 2026
% Reference: 
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
% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addParameter('Alpha',0.05);
p.parse(varargin{:});

alpha = p.Results.Alpha;

freqs = freqs(:)';
F = length(freqs);

[d1,d2] = size(X);

% Detect orientation
if d2 == F
    % epochs × frequencies
elseif d1 == F
    X = X';
else
    error('Frequency vector length does not match any dimension of X.');
end

[E,F] = size(X);

if E < 3
    error('Need at least 3 epochs.');
end

% -----------------------
% Statistics
% -----------------------
mu_mean = mean(X,1);
s       = std(X,0,1);
se      = s ./ sqrt(E);

df = E - 1;

tcrit = tinv(1 - alpha/2, df);

ci_low  = mu_mean - tcrit .* se;
ci_high = mu_mean + tcrit .* se;

mu_ci   = [ci_low' ci_high'];

end
