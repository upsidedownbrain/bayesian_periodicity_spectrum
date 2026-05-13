function [postOdds] = posterior_odds_temperature(bf_psi, varargin)
% This function estimates posterior odds from Bayes factors using an
% empirical prior (population mean) adjusted with a temperature parameter.
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
% bf_psi: [F x M] (frecuencies x signals) in log-scale f: [F x 1] in Hz
% Optional:
p = inputParser;
p.addParameter('EpsClip', 1e-12);   % avoids -Inf/NaN
p.addParameter('Temp', 1);          % prior temperature (default tau=1)
p.parse(varargin{:});
epsc  = p.Results.EpsClip;
tau   = p.Results.Temp;

[F,M] = size(bf_psi);
bf_psi = max(bf_psi, log10(epsc));   % avoids -Inf
postOdds = zeros(F,M);

% --- leave-one-out prior (log-mean-exp) ---
% LOO mean of lo gs = (sum_all - col_i)/(M-1)
sumLogs = sum(bf_psi, 2);              % [F x 1]
for i = 1:M
    log10Prior_i = (sumLogs - bf_psi(:,i)) / max(M-1,1);  % [F x 1]
    % temperature (opcional): prior^{tau} -> tau*logPrior
    log10Prior_i = tau * log10Prior_i;
    postOdds(:,i) = log10Prior_i + bf_psi(:,i);
end
end
