function bf10 = spectra_comparison(spectra_1,spectra_2, tail)
% This function performs a paired t-test using the values of each
% frequency component from two populations of spectra.
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
    
%% Input arguments:
% spectra_1 and spectra_2: m x n matrices, being m each freqeuncy
% component and n each spectrum.
% tail: 'both','right', or 'left' for two or one-tailed tests [both].
% Note that 'right' means spectra_1 > spectra_2.

%% Output arguments:
% bf10: Bayes factor

for j = 1:size(spectra_1,1)
    [bf10(j), ~, ~] = ttest_bf(spectra_1(j,:), spectra_2(j,:),'tail', tail);
end

end
