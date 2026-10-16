function [mu_mean, mu_ci, df, rho1] = ...
    credible_interval_individual_spectrum(X, freqs, varargin)
% CREDIBLE_INTERVAL_INDIVIDUAL_SPECTRUM
%
% Pointwise Bayesian credible intervals for the mean spectrum across
% segments, accounting for first-order serial dependence using an
% AR(1) covariance model.
%
% MODEL
%
%   x ~ N(mu * 1, sigma^2 * R)
%
% where
%
%   R(i,j) = rho1^|i-j|
%
% and the reference prior
%
%   p(mu,sigma^2) proportional to 1/sigma^2
%
% is used.
%
% Under this model, conditional on the estimated AR(1) correlation,
% the marginal posterior of mu is Student-t.
%
% INPUT
%   X     : [segments x frequencies] spectral values
%   freqs : frequency vector
%
% OPTIONS
%   'Alpha' : posterior tail probability (default = 0.05)
%   'SampleDimension' : 1 for segments x frequencies, 2 for
%                       frequencies x segments. Default infers orientation
%                       only when unambiguous; square matrices require it.
%
% OUTPUT
%   mu_mean : [1 x F] posterior location / GLS estimate
%   mu_ci   : [F x 2] pointwise credible intervals
%   df      : posterior degrees of freedom
%   rho1    : estimated lag-1 correlation at each frequency
%
% Authors: Jesús Pardo-Valencia, Guglielmo Foffani
% Version: 2.0.0 
% Date: October 2026
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
% NOTE:
% rho1 is estimated from the observed segment-level data and then treated
% as fixed when computing the posterior for mu and sigma^2. Therefore,
% the resulting credible interval is a plug-in approximation and does not
% integrate uncertainty in the AR(1) correlation parameter.
% Constant frequency columns are rejected: the stated reference prior
% gives an improper posterior when the residual sum of squares is zero.


%% Parse inputs

p = inputParser;

p.addParameter( ...
    'Alpha',0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);

p.addParameter('SampleDimension',[], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ismember(x,[1 2])));

p.parse(varargin{:});

alpha = p.Results.Alpha;


%% Orientation

validateattributes(X,{'numeric'},{'real','finite','nonempty','2d'},mfilename,'X');
validateattributes(freqs,{'numeric'},{'real','finite','nonempty','vector'},mfilename,'freqs');
freqs = freqs(:)';
F = numel(freqs);
sampleDim = p.Results.SampleDimension;
if isempty(sampleDim)
    matches = [size(X,2)==F, size(X,1)==F];
    if all(matches)
        error('credible_interval_individual_spectrum:AmbiguousOrientation', ...
            'Square input is ambiguous. Specify SampleDimension=1 or 2 explicitly.');
    elseif matches(1)
        sampleDim = 1;
    elseif matches(2)
        sampleDim = 2;
    else
        error('credible_interval_individual_spectrum:FrequencySize', ...
            'Frequency vector length does not match the frequency dimension.');
    end
end
if size(X,3-sampleDim) ~= F
    error('credible_interval_individual_spectrum:FrequencySize', ...
        'Frequency vector length does not match the specified frequency dimension.');
end
if sampleDim == 2, X = X.'; end

[N,F] = size(X);


if N < 3
    error('At least 3 segments are required.');
end


%% Preallocate

mu_mean = nan(1,F);
mu_ci   = nan(F,2);

rho1 = nan(1,F);

df = N - 1;


%% Constant vectors

one = ones(N,1);


%% Frequency-wise inference

for f = 1:F

    x = X(:,f);

    if any(~isfinite(x))
        error('X contains non-finite values.');
    end
    if all(x == x(1))
        error('credible_interval_individual_spectrum:DegenerateData', ...
            'Frequency column %d is constant; the reference-prior posterior is improper.',f);
    end


    %% ---------------------------------------------------------------
    % Estimate first-lag autocorrelation
    % ---------------------------------------------------------------

    x0 = x - mean(x);

    if std(x0) == 0

        r1 = 0;

    else

        [r,lags] = xcorr(x0,1,'coeff');

        r1 = r(lags == 1);

        if isempty(r1) || ~isfinite(r1)
            r1 = 0;
        end

    end


    % Avoid nearly singular AR(1) covariance matrices
    r1 = max(min(r1,0.99),-0.99);

    rho1(f) = r1;


    %% ---------------------------------------------------------------
    % AR(1) correlation matrix
    % ---------------------------------------------------------------

    idx = (1:N)';

    R = r1 .^ abs(idx - idx');


    %% ---------------------------------------------------------------
    % Generalized least-squares estimate of mu
    % ---------------------------------------------------------------

    % Prefer matrix division rather than explicitly calculating inv(R)
    Rinv_one = R \ one;
    Rinv_x   = R \ x;

    information_mu = one' * Rinv_one;

    mu_hat = (one' * Rinv_x) / information_mu;


    %% ---------------------------------------------------------------
    % Residual variance under AR(1) model
    % ---------------------------------------------------------------

    residual = x - mu_hat;

    SSE = residual' * (R \ residual);
    if ~isfinite(SSE) || SSE <= 0
        error('credible_interval_individual_spectrum:DegenerateData', ...
            'Frequency column %d has no finite positive residual sum of squares.',f);
    end

    sigma2_hat = SSE / df;


    %% ---------------------------------------------------------------
    % Posterior scale of mu
    %
    % Under p(mu,sigma^2) proportional to 1/sigma^2:
    %
    % mu | X ~ Student-t_df(mu_hat, posterior_scale)
    % ---------------------------------------------------------------

    posterior_scale = ...
        sqrt(sigma2_hat / information_mu);


    %% ---------------------------------------------------------------
    % Pointwise credible interval
    % ---------------------------------------------------------------

    tcrit = tinv(1-alpha/2,df);

    ci_low  = mu_hat - tcrit * posterior_scale;
    ci_high = mu_hat + tcrit * posterior_scale;


    mu_mean(f) = mu_hat;

    mu_ci(f,:) = [ci_low ci_high];

end

end
