% Demonstration of how to simulate the movement of a ship on predefined 
% waves. Assume ship points north.
% Authors: Matheus Bernat & Ludvig Granström 2021
% Last updated: July 2021
%
% Demo 1-3: Sea state 6, long-crested waves, no ship speed,
% Demo 4:   Sea state 3, long-crested waves, no ship speed
% Demo 5:   Sea state 3, long-crested waves, ship speed of 15 knots
% Demo 6:   Sea state 3, short-crested waves, ship speed of 15 knots
% Demo 7:   No waves, stabilization test
%
% To plot and visualize simulations again, run:
% w = load(waveFile).wavesStruct;
% plotShipStates(states, w.tVec, w.Ts, w.beta);
% visualizeSimulation(states, w.waves, w.xVec, w.yVec, w.tVec, face, vert, cogVec);

% Add paths to wave-files, help-files
wavesPath = [pwd, '\wave-files'];
helpFilesPath = [pwd, '\help-files'];
boatPath = [pwd, '\boat-files'];
addpath(wavesPath);
addpath(helpFilesPath);
addpath(boatPath);

%% ------- Create ship with properties below
% Ship:                       HMS Norfolk (hull only)
% Mass:                       2.5e6
% Dimensions:                 137x15x16 
% Vertices initial position:  [30 35 5.55]
% Reference velocity along u: 0 [m/s] (should always be 0, the ship speed is set in the wave file)
% Reference yaw:              0 degrees
% CoG offset:                 [-4.77 0.022 -2]
% Position of heliPad         [40 35 6]
shipStruct.file        = 'filteredNorfolkNew.stl';
shipStruct.M           = 2.5e6;
shipStruct.len         = 137;
shipStruct.width       = 15;
shipStruct.height      = 16;
shipStruct.verticesPos = [30 35 5.55];
shipStruct.refSpeedU   = 0;
shipStruct.refYaw      = 0;
shipStruct.cogOffset   = [-4.77 0.022 -2];
shipStruct.helipadPos  = [40 35 6];
%               [v_u v_v v_w phi th psi w_phi w_th w_psi]'
shipStruct.x0 = [0    0   0  0   0  0   0     0    0]; % Initial state values


%% --------------- Demo #1: Sea state 6 pitch test (waves from the north)
% ------- Use wave with properties below
% Sea state:        6
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       180 degrees (from the front)
% Grid:             300x100
% Time:             0:0.2:500
% Ship speed        0
isPlot = true;
isVisual = true;
waveFile = 'waves__seaState_6__long__beta_3.14__grid_300x100__time_0_0.2_500__U_0.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, isPlot, isVisual);

%% --------------- Demo #2: Sea state 6 roll test (waves from the west)
% ------- Use wave with properties below
% Sea state:        6
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       90 degrees (from the left)
% Grid:             300x100
% Time:             0:0.2:200
% Ship speed        0

waveFile = 'waves__seaState_6__long__beta_1.57__grid_400x200__time_0_0.2_100__U_0.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, true, true);

%% --------------- Demo #3: Sea state 6 pitch & roll test (waves from northwest)
% ------- Use wave with properties below
% Sea state:        6
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       pi degrees (from the left)
% Grid:             300x100
% Time:             0:0.2:500
% Ship speed        0

waveFile = 'waves__seaState_6__long__beta_2.36__grid_300x100__time_0_0.2_200__U_0.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, true, true);

%% --------------- Demo #4: Sea state 3, waves from the north, no speed
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       180 degrees (from the front)
% Grid:             300x100
% Time:             0:0.2:500
% Ship speed        0

waveFile = 'waves__seaState_3__long__beta_3.14__grid_300x100__time_0_0.2_500__U_0.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, true, true);

%% --------------- Demo #5: Sea state 3, waves from the north, 15 knots
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       180 degrees (from the front)
% Grid:             300x100
% Time:             0:0.2:500
% Ship speed        0

waveFile = 'waves__seaState_3__long__beta_3.14__grid_300x100__time_0_0.2_500__U_7.72.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, true, true);

%% --------------- Demo #6: Sea state 3, waves from the north, 15 knots
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): short-crested (multi-directional)
% Wave angle:       180 degrees (from the front)
% Grid:             300x100
% Time:             0:0.2:500
% Ship speed        0

waveFile = 'waves__seaState_3__short__beta_3.14__grid_300x100__time_0_0.2_500__U_7.72.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, true, true);

%% --------------- Demo #7: Stabilization test
% ------- Use wave with properties below
% Sea state:        3
% Wave type :       No waves
% Grid:             300x100
% Time:             0:0.2:50
% Ship speed        0

waveFile = 'waves_seaState_3_long_beta_3.14_grid_300x100_time_0_0.2_500.mat';
[states, face, vert, cogVec] = simulateShip(waveFile, shipStruct, true, true);




