function results = periodicity_analysis(s, Fs, window_BF, method, psi_non_osc, opts)
%% Input arguments
% s: Input signal (column vector)
% Fs: Sampling frequency
% window_BF: Width of the segments (in samples) in which the input
% signal is going to be divided for posterior statistics.
% method: Short-time Fourier transform ('STFT') or Continuous wavelet
% transform ('CWT').
% psi_non_osc: Structure file of asymptotic periodicity spectrum of
% non-oscillatory activity 𝛹𝑛𝑜𝑛−𝑜𝑠𝑐𝑖𝑙𝑙𝑎𝑡𝑜𝑟𝑦(𝑓)
%
% Optional:
% window: For STFT. Length (in samples) of the Hann window.  
% overlapLength: For STFT. Overlapping (in samples) of the Hann window.
% FFTLength: For STFT. Length (in samples) of FFT.
% VoicesPerOctave: For CWT. Voices per octave.
% CI: Credible intervals for the normalized periodicity spectrum.
% BF_corr: Correction of effective sample size by the first
% autocorrelation lag for Bayesian periodicity spectrum.
%
% Example:
% results = periodicity_analysis(signal, 1000, 2000, "STFT", ...
% overlapLength=0.5, FFTLength=2048);

%% Output arguments
% results: struct type

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

arguments
    s
    Fs
    window_BF
    method {mustBeMember(method, ["STFT", "CWT"])}
    psi_non_osc
    opts.window (1,1) double {mustBePositive, mustBeInteger} = Fs
    opts.overlapLength (1,1) double {mustBeNonnegative} = Fs-1
    opts.FFTLength (1,1) double {mustBePositive, mustBeInteger} = Fs
    opts.VoicesPerOctave (1,1) double {mustBePositive, mustBeInteger} = 16
    opts.CI (1,1) logical = false
    opts.BF_corr (1,1) logical = true
end

%% Variables in columns
if size(s,1) == 1
    s = s';
end

%% 1. Frequency decomposition
% Optional: use analytic signal to avoid negative frequencies mirroring
% s = hilbert(s);
if strcmp(method, 'STFT')
    [tfr, freqs, ~] = stft(s, Fs, 'Window', hann(opts.window, 'periodic'), 'OverlapLength',...
        opts.overlapLength, 'FFTLength', opts.FFTLength);
    % Keep only positive frequencies
    pos = freqs > 0 & freqs < (Fs/2);
    tfr = tfr(pos, :);
    freqs   = freqs(pos);
elseif strcmp(method, 'CWT')
    [tfr, freqs] = cwt(s, 'amor', Fs, 'voicesPerOctave', opts.VoicesPerOctave);
    % Frequencies above 1 Hz
    tfr(freqs<1,:) = [];
    freqs(freqs<1) = [];
    FREQ = centfrq('morl');
    scale = []; scale = (FREQ./freqs);
else
    error('Not valid. Use STFT or CWT.');
end
results.freqs = freqs;

    %% 2. Instantaneous frequency (phase-vocoder correction)
if strcmp(method, 'STFT')
    % Hop in samples
    hop = opts.window - opts.overlapLength;

    % F in Hz -> angular frequency per sample of each bin
    omega_k = 2*pi*freqs/Fs;

    % Phase per bin and frame
    phi = angle(tfr);

    % Phase increment between frames
    dphi = diff(phi, 1, 2);

    % Expected phase advance per bin between frames
    expected = omega_k * hop;

    % Residual and projection at the main interval (-pi, pi]
    resid = dphi - expected;                      
    dphi_corr = mod(resid + pi, 2*pi) - pi;

    % Instantaneous frequency per bin
    omega_inst = omega_k + dphi_corr / hop;
    f_inst = omega_inst * (Fs/(2*pi));

elseif strcmp(method, 'CWT')
    phi  = angle(tfr);                 % [Nfreq x Ntime]
    dphi = diff(phi,1,2);              % hop=1

    omega0 = 2*pi*freqs/Fs;                % [Nfreq x 1], expected progress (rad/muestra)
    resid  = dphi - omega0;           
    dphi_corr = mod(resid + pi, 2*pi) - pi;   % to (-pi,pi]

    omega_inst = omega0 + dphi_corr;          % rad/sample
    f_inst = (Fs/(2*pi)) * omega_inst;  % Hz
end

f_inst = f_inst';

%% 3. Power spectrum
segments = floor(size(f_inst,1)/(window_BF));
for z = 1:1:segments
    if strcmp(method, 'STFT')
        power(:,z) = mean(abs(tfr(:,(z-1)*(window_BF)+1:z*((window_BF)))).^2, 2)'; % Measured as the inverse logarithm of the variance of the instantaneous frequency
    elseif strcmp(method, 'CWT')
        power(:,z) = mean(abs(tfr(:,(z-1)*(window_BF)+1:z*((window_BF)))).^2 .* scale, 2)'; % % Measured as the inverse logarithm of the variance of the instantaneous frequency
    end
end
results.power_sp = mean(power,2);

%% 4. Normalized periodicity spectrum
if strcmp(method, 'STFT')
    for z = 1:segments
        psi(:,z) = 1 ./ std(f_inst((z-1)*(window_BF)+1:z*(window_BF),:)) ./ psi_non_osc.mean_psi_non_osc_fitted'; % Measured as the inverse logarithm of the variance of the instantaneous frequency
    end
    results.psi_norm = mean(psi, 2);
elseif strcmp(method, 'CWT')
    for z = 1:segments
        psi(:,z) = 1./std(f_inst((z-1)*(window_BF)+1:z*(window_BF),:)) ./ psi_non_osc.mean_psi_non_osc_fitted'; % Measured as the inverse logarithm of the variance of the instantaneous frequency
    end
    results.psi_norm = mean(psi, 2);
end

%% 5. Credible intervals for normalized periodicity spectrum (optional)
if opts.CI
    [mu_mean, mu_ci, df] = credible_interval_individual_spectrum(psi, freqs);
    results.mu_mean = mu_mean;
    results.mu_ci   = mu_ci;
    results.df      = df;
end

%% 6. Bayesian periodicity spectrum
psi = psi - 1;
if opts.BF_corr
    R_value = [];
    for z = 1:1:length(freqs)
        [r, ~] = xcorr(psi(z,:)-mean(psi(z,:)),'normalized');
        R_value(z) = r(segments + 1);

        bf_psi(z) = ttest_bf(psi(z,:), 'rho1', R_value(z), 'tail', 'right');
        bf_psi(z) = log10(bf_psi(z));
    end
else
    for z = 1:1:length(freqs)
        bf_psi(z) = ttest_bf(psi(z,:), 'tail', 'right');
        bf_psi(z) = log10(bf_psi(z));
    end
end
results.bf_psi = bf_psi;

end
