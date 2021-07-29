function Hs=getSignificantWaveHeight(seaState)
% Get Hs given sea state. Reference: Fossen's handbook pg. 204.
    if seaState == 1 
        Hs = (0.3 - 0) * rand;         % Hs between 0 and 0.3
    elseif seaState == 2 
        Hs = 0.3 + (0.6 - 0.3) * rand; % Hs between 0.3 and 0.6
    elseif seaState == 3 
        Hs = 1 + (2 - 1) * rand;       % Hs between 0.6 and 2
    elseif seaState == 4
        Hs = 2 + (3 - 2) * rand;       % Hs between 2 and 3
    elseif seaState == 5
        Hs = 3 + (4 - 3) * rand;       % Hs between 3 and 4
    elseif seaState == 6
        Hs = 4 + (6 - 4) * rand;       % Hs between 4 and 6
    elseif seaState == 7
        Hs = 6 + (9 - 6) * rand;       % Hs between 6 and 9
    elseif seaState == 8
        Hs = 9 + (14 - 9) * rand;      % Hs between 9 and 14
    end
end