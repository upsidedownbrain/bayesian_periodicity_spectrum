function results = periodicity_analysis(s, Fs, window_BF, method, psi_non_osc, opts)
%% PERIODICITY_ANALYSIS
%
% Estimates:
%   1. Time-frequency decomposition
%   2. Instantaneous frequency
%   3. Frequency-normalized power spectrum
%   4. Normalized periodicity spectrum
%   5. Optional Bayesian credible intervals
%   6. Bayesian periodicity spectrum
%
% INPUTS
% s:
%   Input signal.
%
% Fs:
%   Sampling frequency in Hz.
%
% window_BF:
%   Width, in samples, of the non-overlapping segments used to estimate
%   periodicity and Bayesian periodicity.
%   STFT boundaries are mapped to the next available frame when window_BF
%   is not a multiple of the hop. Each segment needs at least two frames.
%
% method:
%   "STFT" or "CWT".
%
% psi_non_osc:
%   Structure containing:
%
%       psi_non_osc.mean_psi_non_osc_fitted
%
%   with the expected periodicity of non-oscillatory activity at each
%   analysed frequency.
%
% OPTIONAL PARAMETERS
% window:
%   STFT Hann-window length in samples.
%
% overlapLength:
%   STFT overlap in samples.
%
% FFTLength:
%   STFT FFT length.
%
% VoicesPerOctave:
%   Number of voices per octave for CWT.
%
% CI:
%   Compute pointwise Bayesian credible intervals for the normalized
%   periodicity spectrum using an AR(1) covariance model and the
%   reference prior p(mu,sigma^2) proportional to 1/sigma^2.
%
% BF_corr:
%   Apply AR(1)-based effective-sample-size correction to the Bayesian
%   periodicity spectrum.
%
% OUTPUT
% results:
%   Structure containing the estimated spectra and diagnostic variables.
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


arguments
    s
    Fs
    window_BF
    method {mustBeMember(method, ["STFT", "CWT"])}
    psi_non_osc

    opts.window (1,1) double ...
        {mustBePositive, mustBeInteger} = Fs

    opts.overlapLength (1,1) double ...
        {mustBeNonnegative, mustBeInteger} = Fs-1

    opts.FFTLength (1,1) double ...
        {mustBePositive, mustBeInteger} = Fs

    opts.VoicesPerOctave (1,1) double ...
        {mustBePositive, mustBeInteger} = 16

    opts.CI (1,1) logical = false

    opts.BF_corr (1,1) logical = true
end


%% ========================================================================
% 0. Input preparation
% ========================================================================

s = s(:);

if window_BF <= 0 || mod(window_BF,1) ~= 0
    error('window_BF must be a positive integer.');
end


%% ========================================================================
% 1. Frequency decomposition
% ========================================================================

switch method

    case "STFT"

        hop = opts.window - opts.overlapLength;

        if hop <= 0
            error('STFT overlapLength must be smaller than window length.');
        end

        [tfr, freqs, ~] = stft( ...
            s, ...
            Fs, ...
            'Window', hann(opts.window,'periodic'), ...
            'OverlapLength', opts.overlapLength, ...
            'FFTLength', opts.FFTLength);

        keep = freqs > 0 & freqs < Fs/2;

        tfr   = tfr(keep,:);
        freqs = freqs(keep);


    case "CWT"

        [tfr, freqs] = cwt( ...
            s, ...
            'amor', ...
            Fs, ...
            'VoicesPerOctave', opts.VoicesPerOctave);

        keep = freqs >= 1;

        tfr   = tfr(keep,:);
        freqs = freqs(keep);

end


freqs = freqs(:);

results.freqs = freqs;
results.tfr   = tfr;


%% ========================================================================
% 2. Instantaneous frequency
% ========================================================================

switch method

    case "STFT"

        hop = opts.window - opts.overlapLength;

        omega_k = 2*pi*freqs/Fs;

        phi = angle(tfr);

        dphi = diff(phi,1,2);

        expected = omega_k * hop;

        resid = dphi - expected;

        dphi_corr = mod(resid + pi,2*pi) - pi;

        omega_inst = omega_k + dphi_corr/hop;

        f_inst = omega_inst * Fs/(2*pi);


    case "CWT"

        phi = angle(tfr);

        dphi = diff(phi,1,2);

        omega0 = 2*pi*freqs/Fs;

        resid = dphi - omega0;

        dphi_corr = mod(resid + pi,2*pi) - pi;

        omega_inst = omega0 + dphi_corr;

        f_inst = omega_inst * Fs/(2*pi);

end


% time x frequency
f_inst = f_inst.';

results.f_inst = f_inst;


%% ========================================================================
% 3. Segment definition
% ========================================================================

nTime = size(f_inst,1);

sampleHop = 1;
if method == "STFT", sampleHop = hop; end
segments = floor(nTime*sampleHop/window_BF);

if segments < 3
    error([ ...
        'At least three complete segments are required. ' ...
        'Decrease window_BF or provide a longer signal.' ...
        ]);
end


segmentEdges = ceil((0:segments)*window_BF/sampleHop);
if any(diff(segmentEdges) < 2)
    error('periodicity_analysis:InsufficientSegmentFrames', ...
        'Each segment requires at least two instantaneous-frequency frames. Increase window_BF or decrease the STFT hop.');
end
nUsed = segmentEdges(end);

f_inst_used = f_inst(1:nUsed,:);


results.n_segments         = segments;
results.n_timepoints_used  = nUsed;


%% ========================================================================
% 4. Power spectrum
% ========================================================================

nFreq = numel(freqs);

power_seg = nan(nFreq,segments);


% Frequency normalization for CWT power
%
% White-noise simulations showed that multiplication by 1/f removes
% the systematic frequency dependence of mean squared CWT magnitude.
%
% This quantity should therefore be interpreted as frequency-normalized
% CWT power rather than as an absolutely calibrated PSD.

if method == "CWT"
    cwt_scale = 1 ./ freqs;
end


for z = 1:segments

    idx = (segmentEdges(z)+1):segmentEdges(z+1);


    switch method

        case "STFT"

            % Mean squared magnitude:
            % E[|X|^2]
            power_seg(:,z) = ...
                mean(abs(tfr(:,idx)).^2,2);


        case "CWT"

            % Mean squared CWT magnitude with frequency normalization
            power_seg(:,z) = ...
                mean(abs(tfr(:,idx)).^2,2) .* cwt_scale;

    end

end


results.power_segments = power_seg;
results.power_sp       = mean(power_seg,2);


%% ========================================================================
% 5. Normalized periodicity spectrum
% ========================================================================

psi = nan(nFreq,segments);

psi0 = psi_non_osc.mean_psi_non_osc_fitted(:);


if numel(psi0) ~= nFreq
    error([ ...
        'psi_non_osc.mean_psi_non_osc_fitted must contain one value ' ...
        'for every analysed frequency.' ...
        ]);
end


for z = 1:segments

    idx = (segmentEdges(z)+1):segmentEdges(z+1);

    % Standard deviation of instantaneous frequency within segment
    sigma_IF = std( ...
        f_inst_used(idx,:), ...
        0, ...
        1).';

    % Raw periodicity
    psi_raw = 1 ./ sigma_IF;

    % Normalized periodicity
    psi(:,z) = psi_raw ./ psi0;

end


results.psi_segments = psi;
results.psi_norm     = mean(psi,2);


%% ========================================================================
% 6. Bayesian credible intervals
% ========================================================================
%
% Pointwise credible intervals for the mean normalized periodicity
% across segments are obtained using:
%
%   x ~ N(mu * 1, sigma^2 * R)
%
% where R is an AR(1) correlation matrix:
%
%   R(i,j) = rho1^|i-j|
%
% and the reference prior:
%
%   p(mu,sigma^2) proportional to 1/sigma^2.
%
% rho1 is estimated independently at each frequency.
%

if opts.CI

    [ ...
        mu_mean, ...
        mu_ci, ...
        df_ci, ...
        rho1_ci ...
    ] = credible_interval_individual_spectrum( ...
            psi, ...
            freqs, ...
            'SampleDimension',2);


    results.mu_mean = mu_mean;
    results.mu_ci   = mu_ci;

    % Degrees of freedom of the posterior Student-t
    results.df_ci = df_ci;
    results.df = df_ci;

    % Estimated AR(1) correlation used for each frequency
    results.rho1_ci = rho1_ci;

end


%% ========================================================================
% 7. Bayesian periodicity spectrum
% ========================================================================
%
% H0:
%
%       Psi_norm = 1
%
% H+:
%
%       Psi_norm > 1
%
% Therefore subtract 1 so that the point null corresponds to zero.
%

psi_test = psi - 1;


bf10_values = nan(nFreq,1);
log10_bf    = nan(nFreq,1);

rho1_values = nan(nFreq,1);

T_values = nan(nFreq,1);
N_values = nan(nFreq,1);
p_values = nan(nFreq,1);


for z = 1:nFreq

    % Segment-level normalized periodicity values
    x = psi_test(z,:).';


    %% --------------------------------------------------------------------
    % Effective-sample-size correction
    % ---------------------------------------------------------------------

    if opts.BF_corr

        x0 = x - mean(x);


        if std(x0) == 0

            rho1 = 0;

        else

            [r,lags] = xcorr( ...
                x0, ...
                1, ...
                'coeff');

            rho1 = r(lags == 1);


            if isempty(rho1) || ~isfinite(rho1)
                rho1 = 0;
            end

        end


        rho1_values(z) = rho1;


        %% ---------------------------------------------------------------
        % Directional Bayesian t-test
        % ---------------------------------------------------------------
        %
        % ttest_bf applies:
        %
        %   N_eff = N * (1-rho1)/(1+rho1)
        %
        % and:
        %
        %   T_eff = T * sqrt(N_eff/N)
        %
        % and evaluates the directional BF using:
        %
        %   - noncentral-t likelihood
        %   - normalized positive half-Cauchy prior
        %   - Cauchy scale = 1/sqrt(2)
        %   - H+ : delta > 0
        %

        [ ...
            bf10_values(z), ...
            p_values(z), ...
            T_values(z), ...
            N_values(z), ...
            log10_bf(z) ...
        ] = ttest_bf( ...
                x, ...
                'rho1',rho1, ...
                'tail','right');


    else

        rho1_values(z) = 0;


        [ ...
            bf10_values(z), ...
            p_values(z), ...
            T_values(z), ...
            N_values(z), ...
            log10_bf(z) ...
        ] = ttest_bf( ...
                x, ...
                'tail','right');

    end


    % ttest_bf returns log10 BF directly, preserving finite evidence
    % even if the raw BF exceeds floating-point range.

end


%% ========================================================================
% 8. Store Bayesian outputs
% ========================================================================

% Original toolbox output:
% log10 directional Bayes factor
results.bf_psi = log10_bf.';


% Raw BF+0
results.bf10_psi = bf10_values;


% Diagnostics
results.rho1     = rho1_values;
results.T_bf     = T_values;
results.N_eff_bf = N_values;
results.p_bf     = p_values;


end
