function waves=simulateWaves(seaState, xVec, yVec, beta, tVec, U, lambda, muVec, dmu)
% SIMULATEWAVES(seaState, xVec, yVec, beta, tVec, U , lambda, muVec, dmu) 
% Takes in sea state and plots a wave height. Uses the Bretschneider spectrum.
%
% Inputs:
%   - seaState: integer in interval [1, 9].
%   - xVec:     1xN vector. E.g. xVec = linspace(-50, 50, 100).
%   - yVec:     1xN vector. E.g. yVec = linspace(-50, 50, 100).
%   - beta:     direction of main wave in rad.
%   - tVec:     time vector. E.g.: tVec = 0:0.1:50.
%   - U:        speed of the ship [m/s]
%   - muVec:    angle directions vector. E.g.: -pi/2:dmu:pi/2. If muVec=[],
%               then the waves generated will be unidirectional (long-
%               crested).
%   - dmu:      direction interval taken in muVec. Used if muVec ~= [].
%   - lambda:   waveLength. If a sea with infinite depth is assumed, set
%               lambda to [].
%
% Ouput:
%   - waves:    if ~is3d, then size(waves) = (1, length(tVec)), where waves
%               is equal to the wave height in some point in the sea. If  
%               is3d, then size(waves) = (length(xVec), length(yVec), 
%               length(tVec)), where waves will then represent the wave
%               height of all points (x, y) over a 2D grid defined by xVec
%               and yVec over an interval of time tVec.
tic; 
disp('Creating waves...');
g = 9.81;

% Get significant wave height
Hs = getSignificantWaveHeight(seaState);

% Create Bretschneider spectrum given Hs and plot spectrum
dw = 0.1;
wVec = (dw/2:dw:3)';

A = 8.1 * 1e-3 * g^2; % constant, eq 8.54
B = 3.11 / (Hs^2);    % eq 8.55
specType = 1; % Bretschneider (@ Fossen pg 203), from kravspec 
S = wavespec(specType, [A, B], wVec, 0);
S(1) = 0; % the first element is NaN for some reason

waves = zeros(length(yVec), length(xVec), length(tVec));

% Get the set of frequencies, directions and phases that will be used to
% generate waves for all points in the grid:
waveFrequencies = zeros(1, length(wVec));
for k=1:length(waveFrequencies)        
    waveFrequencies(k) = wVec(k) - dw/2 + dw * rand;
end

if nargin > 7 % Short-crested wave
    waveDirections = zeros(1, length(muVec));
    for i=1:length(waveDirections)
        waveDirections(i) = muVec(i) - dmu/2 + dmu * rand;            
    end
    sizeWaveDirections = length(waveDirections);
else
    sizeWaveDirections = 1;
end

wavePhases = zeros(sizeWaveDirections, length(waveFrequencies));
for k=1:length(waveFrequencies) 
    for i=1:sizeWaveDirections
        wavePhases(k, i) = 2 * pi * rand;
    end 
end

% Get the wave heights for all coordinates (x,y) for all times in tVec
for yIdx = 1:length(yVec)
    y = yVec(yIdx);
    for xIdx=1:length(xVec)
        x = xVec(xIdx);

        % For each coordinate (x,y), sum the contribution of all the waves
        % to get the final wave amplitude at that point = sumOfWaves.
        sumOfWaves = 0;

        for k=1:length(waveFrequencies)
            w_k = waveFrequencies(k);
            
            % Long-crested
            if nargin < 8
                if nargin > 6
                    % Lambda is passed as an argument
                    coeff = 2 * pi / lambda; 
                else
                    % Infinite depth sea assumed
                    coeff = w_k ^ 2 / g;
                end
                e_k = wavePhases(k);
                amp = sqrt(2 * S(k) * dw);
                wave = amp * cos(coeff * ((x+U*tVec) * cos(-beta) ...
                               + y * sin(-beta))  ...
                                          - w_k * tVec + e_k);
                sumOfWaves = sumOfWaves + wave;
                
            % Short-crested
            else              
                for i=1:length(waveDirections)
                    e_ik = wavePhases(1, k);
                    mu_i = waveDirections(i);

                    amp = sqrt(2 * S(k) * spread(mu_i) * dw * dmu);
                    wave = amp * cos(w_k^2/g * ((x+U*tVec) * cos(mu_i - beta) ...
                                   + y * sin(mu_i - beta))  ...
                                              - w_k * tVec + e_ik);
                    sumOfWaves = sumOfWaves + wave;                
                end 
            end
        end
        waves(yIdx, xIdx, :) = sumOfWaves;
    end
end   
disp('Done creating waves!');
toc;


% ----------------------- Help functions ----------------------------------

function spreadFunction=spread(mu)
    if (mu >= -pi/2 && mu <= pi/2)
        spreadFunction = (2 / pi) * cos(mu)^2;
    else
        spreadFunction = 0;
    end
end
end
