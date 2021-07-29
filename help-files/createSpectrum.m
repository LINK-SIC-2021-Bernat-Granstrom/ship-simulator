function S = createSpectrum(seaState,wVec) 
% CREATESPECTRUM creates the Bretschneider spectrum given sea state and
% a frequency vector 
% Inputs:
%   - seaState: Integer in interval [1, 9].
%   - wVec:     Vector containing the frequencies. e.g (0:0.01:3)'
%Ouput:
%   - S:      Vector containing the spectral density for each frequency 
g = 9.81;

% Get significant wave height
Hs = getSignificantWaveHeight(seaState);

% Create Bretschneider spectrum
A = 8.1 * 1e-3 * g^2; % constant, eq 8.54
B = 3.11 / (Hs^2);    % eq 8.55
specType = 1; % Bretschneider (@ Fossen pg 203), from kravspec 
S = wavespec(specType, [A, B], wVec, 0);
S(1) = 0; % the first element is NaN for some reason

end