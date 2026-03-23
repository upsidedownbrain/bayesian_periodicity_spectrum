function results = non_oscillatory_periodicity_spectrum(signal_duration, Fs, window_BF, method, opts)
% This function obtains the asymptotic periodicity spectrum of non-oscillatory activity 𝛹𝑛𝑜𝑛−𝑜𝑠𝑐𝑖𝑙𝑙𝑎𝑡𝑜𝑟𝑦(𝑓)

%% Input arguments
% signal_duration: Length (in samples). If there is more than one
% signal, signal_duration is an array containing the minimum and maximum length of the signals of interest.
% Fs: Sampling frequency.
% window_BF: Length of the segments (in samples) in which the input
% signal is going to be divided for posterior statistics (Bayesian periodicity spectrum).
% method: Short-time Fourier transform ('STFT') or Continuous wavelet
% transform ('CWT').

% Optional:
% window: For STFT. Length (in samples) of the Hann window.  
% overlapLength: For STFT. Overlapping (in samples) of the Hann window.
% FFTLength: For STFT. Length (in samples) of FFT.
% VoicesPerOctave: For CWT. Voices per octave.
% N_simulations: Number of simulations to obtain periodicity of
% non-oscillatory activity

% Example:
% [psi, f] = non_oscillatory_periodicity_spectrum([3000, 4000], 1000, 'STFT', ...
% window=3, ...
% overlapLength=1.5, ...
% FFTLength=4096, ...
% VoicesPerOctave=48);

%% Output arguments
% Results: struct type

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
    signal_duration
    Fs
    window_BF
    method {mustBeMember(method, ["STFT", "CWT"])}
    opts.window (1,1) double {mustBePositive, mustBeInteger} = round(1*Fs)
    opts.hop (1,1) double {mustBePositive, mustBeInteger} = 1
    opts.FFTLength (1,1) double {mustBePositive, mustBeInteger} = round(1*Fs)
    opts.VoicesPerOctave (1,1) double {mustBePositive, mustBeInteger} = 16
    opts.N_simulations (1,1) double {mustBePositive, mustBeInteger} = 10000
end

rng('shuffle');

if strcmp(method, 'STFT')
elseif strcmp(method, 'CWT')
else
    error('Not valid. Use STFT or CWT.');
end
% We calculate the periodicity of 10000 non-oscillatory signals
if length(signal_duration) == 1
    s_l = round(ones(10,1)*signal_duration);
elseif length(signal_duration) == 2
    if signal_duration(1) > signal_duration(2)
        error('The first element must be lower than the second one.')
    end
    s_l = round(signal_duration(1):(signal_duration(2)-signal_duration(1))/9:signal_duration(2));
else
    error('"signal_duration" cannot have more than 2 elements.')
end

non_osc = []; psi_non_osc = [];
for j1 = 1:1:length(s_l)
    cn = dsp.ColoredNoise(0,s_l(j1),round(opts.N_simulations)/length(s_l));
    s = cn();
    s = (s - mean(s)) ./ std(s);
    for i = 1:1:size(s,2)
        %% 1. Frequency decomposition
        if strcmp(method, 'STFT')
            [tfr, freqs, ~] = stft(s(:,i), Fs, 'Window', hann(opts.window, 'periodic'), 'OverlapLength',...
                opts.window - opts.hop, 'FFTLength', opts.FFTLength);
            % Keep only positive frequencies (analytic signal has no negatives)
            pos = freqs > 0 & freqs < (Fs/2);                   % drop DC and Nyquist
            tfr = tfr(pos, :);
            freqs   = freqs(pos);

        elseif strcmp(method, 'CWT')
            [tfr, freqs] = cwt(s(:,i), 'amor', Fs, 'voicesPerOctave', opts.VoicesPerOctave);
            tfr(freqs<1,:) = [];
            freqs(freqs<1) = [];
            FREQ = centfrq('morl');
            scale = []; scale = (FREQ./freqs);

        end
        results.freqs = freqs;

        %% 2. Instantaneous frequency (phase-vocoder correction)
        if strcmp(method, 'STFT')
            % Hop in samples
            hop = opts.hop;

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
            dphi = diff(phi,1,2);              % rad por muestra (hop=1)

            omega0 = 2*pi*freqs/Fs;                % [Nfreq x 1], avance esperado (rad/muestra)
            resid  = dphi - omega0;           
            dphi_corr = mod(resid + pi, 2*pi) - pi;   % a (-pi,pi]

            omega_inst = omega0 + dphi_corr;          % rad/muestra
            f_inst = (Fs/(2*pi)) * omega_inst;  % Hz

        end
        f_inst_final = f_inst';

        %% 3. Periodicity
        segments = floor(size(f_inst_final,1)/window_BF);
        for z = 1:segments
            psi = 1./std(f_inst_final((z-1)*(window_BF)+1:z*(window_BF),:)); % Measured as the inverse logarithm of the variance of the instantaneous frequency
            non_osc(j1).value(:,i,z) = psi;
        end
    end
    if j1 > 1
        lim = length(psi_non_osc(:,1));
        psi_non_osc(:,j1) = mean(mean(non_osc(j1).value(1:lim,:,:),3),2);
    else
        psi_non_osc(:,j1) = mean(mean(non_osc(j1).value,3),2);
        f_psi_non_osc_final = freqs;
    end
end

results.mean_psi_non_osc = mean(psi_non_osc,2);
if strcmp(method, 'STFT')
    results.fit_mean = mean(results.mean_psi_non_osc(2:end-1));
    results.mean_psi_non_osc_fitted = mean(results.mean_psi_non_osc(2:end-1)).*ones(size(results.mean_psi_non_osc));

elseif strcmp(method, 'CWT')
    x = log(results.freqs);
    y = log(results.mean_psi_non_osc);

    results.fit_mean = fit(x, y, 'poly6');
    results.mean_psi_non_osc_fitted = exp(results.fit_mean.p1.*log(freqs).^6 +...
        results.fit_mean.p2.*log(freqs).^5 + results.fit_mean.p3.*log(freqs).^4 +...
        results.fit_mean.p4.*log(freqs).^3 + results.fit_mean.p5.*log(freqs).^2 +...
        results.fit_mean.p6.*log(freqs) + results.fit_mean.p7);

end
end
