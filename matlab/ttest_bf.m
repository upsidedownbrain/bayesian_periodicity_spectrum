function [bf10, pValue, T, N] = ttest_bf(X, varargin)
% TTEST_BF_DIRECTIONAL Bayes Factor (directional) and posterior for one-sample / paired t-tests
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
% Usage:
% [bf10,pValue,T] = ttest_bf_directional(X, ...)
% Optional Parm/Value pairs:
%   'tail'  - 'both','right','left'      (default 'both')
%   'scale' - r of Cauchy prior (default sqrt(2)/2)
%   'T'     - supply T value directly
%   'N'     - sample size(s)
%   'rho1'  - first-lag autocorrelation for ESS
%
% Output:
%   bf10    - Bayes factor H1 (δ>0) / H0
%   pValue  - frequentist p-value
%   T       - t-statistic
%   N       - effective sample size

if isnumeric(X)
    if mod(numel(varargin),2)==0
        Y = 0;            % one-sample vs 0
        parms = varargin;
    else
        if numel(varargin)>1
            parms = varargin(2:end);
        else
            parms = {};
        end
        Y  = varargin{1}; % could be scalar M or vector for paired
    end
else
    % call with 'T','N',...
    parms = cat(2,X,varargin);
    X=[]; Y=[];
end

%% Parse inputs
p = inputParser;
p.addParameter('tail','both',@(x) ischar(x) && ismember(upper(x),{'BOTH','RIGHT','LEFT'}));
p.addParameter('scale',sqrt(2)/2,@isnumeric);
p.addParameter('T',[],@isnumeric);
p.addParameter('N',[],@isnumeric);
p.addParameter('rho1',[],@isnumeric);
p.parse(parms{:});

tail = lower(p.Results.tail);
r    = p.Results.scale;
rho1 = p.Results.rho1;

%% Compute frequentist T if not provided
if isempty(p.Results.T)
    [~,pValue,~,stats] = ttest(X,Y,'tail',tail);
    T  = stats.tstat;
    df = stats.df;
    N  = numel(X);

    % --- ESS correction if rho1 ---
    if ~isempty(rho1)
        r1 = max(min(rho1,0.99),-0.99);
        Neff = max(2, min(N, N*(1-r1)/(1+r1)));
        df   = Neff-1;
        N    = Neff;
    end
else
    T = p.Results.T;
    N = p.Results.N;
    if numel(N)==2
        df = sum(N)-2;
        N  = prod(N)/sum(N);
    else
        df = N-1;
    end
end

% Frequentist p-value
switch tail
    case 'both'
        pValue = 2*(1 - tcdf(abs(T), df));
    case 'right'
        pValue = 1 - tcdf(T, df);
    case 'left'
        pValue = tcdf(T, df);
end

%% --- Bayes factor directional (half-Cauchy prior δ>0) ---
% Marginal likelihood of data given delta
t_likelihood = @(delta) ((1 + (T - sqrt(N)*delta).^2/df).^(-(df+1)/2));

% Half-Cauchy prior density for delta>0
halfCauchy = @(delta) (2/pi/r) * (1 ./ (1 + (delta/r).^2)); % normalized over δ>0

% Integrand for directional BF
integrand = @(delta) t_likelihood(delta) .* halfCauchy(delta);

% Compute marginal likelihood for H1 (δ>0)
ml_H1 = integral(integrand, 0, Inf, 'RelTol',1e-6,'AbsTol',1e-10);

% Likelihood for point null δ=0
ml_H0 = t_likelihood(0);

% Bayes factor H1/H0 (directional)
bf10 = ml_H1 / ml_H0;

end
