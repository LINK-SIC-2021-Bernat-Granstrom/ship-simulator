function Hs=getSignificantWaveHeight(seaState)
% GETSIGNIFICANTWAVEHEIGHT Returns the significant wave height Hs given the
% sea state from 0 to 8. Reference: Fossen's handbook pg. 204.
    if seaState == 0 
        Hs = 0;       
    elseif seaState == 1 
        Hs = (0.3 - 0) * rand;         % Hs between 0 and 0.3 m
    elseif seaState == 2 
        Hs = 0.3 + (0.6 - 0.3) * rand; % Hs between 0.3 and 0.6 m
    elseif seaState == 3 
        Hs = 1 + (2 - 1) * rand;       % Hs between 0.6 and 2 m
    elseif seaState == 4
        Hs = 2 + (3 - 2) * rand;       % Hs between 2 and 3 m
    elseif seaState == 5
        Hs = 3 + (4 - 3) * rand;       % Hs between 3 and 4 m
    elseif seaState == 6
        Hs = 4 + (6 - 4) * rand;       % Hs between 4 and 6 m
    elseif seaState == 7
        Hs = 6 + (9 - 6) * rand;       % Hs between 6 and 9 m
    elseif seaState == 8
        Hs = 9 + (14 - 9) * rand;      % Hs between 9 and 14 m
    end
end