% Demonstration of how to create wave files with different wave properties.
% Created by: Matheus Bernat & Ludvig Granström in 2021

% Add paths to wave-files, help-files
helpFilesPath = [pwd, '\help-files'];
addpath(helpFilesPath);

% To visualize waves, run:
% waveFile = 'waves_seaState_6_long_beta_3.14_grid_300x100_time_0_0.2_500_U_0.mat';
% w = load(waveFile).wavesStruct;
% sea = surf(w.xVec, w.yVec, w.waves(:, :, 1));
% axis([w.xVec(1) w.xVec(end) w.yVec(1) w.yVec(end) -5 5]);
% pause(0.1);
% for t=2:length(w.tVec)-1
%     sea.ZData = w.waves(:, :, t);
%     pause(0.02);
%     drawnow;
% end

%% Wave #1: Sea state 6 wave coming from bow (front)
% ------- Use wave with properties below
% Sea state:        6
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       180 degrees
% Grid:             200x100
% Time:             0:0.2:500
% Relative speed    0 m/s

wavesStruct.seaState = 6;
wavesStruct.beta     = pi;
wavesStruct.xVec     = linspace(0, 299, 300);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:200;
wavesStruct.U        = 0;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec,wavesStruct.U);
wavesStruct.waveType = 'long';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Wave #2: Sea state 6 wave coming from port (left)
% ------- Use wave with properties below
% Sea state:        6
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       pi/2 degrees
% Grid:             300x100
% Time:             0:0.2:100
% Relative speed    0 m/s

wavesStruct.seaState = 6;
wavesStruct.beta     = pi / 2;
wavesStruct.xVec     = linspace(0, 299, 300);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:100;
wavesStruct.U        = 0;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec, wavesStruct.U);
wavesStruct.waveType = 'long';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Wave #3: Sea state 6 wave coming from off the port bow 
% ------- Use wave with properties below
% Sea state:        6
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       3*pi/4 degrees
% Grid:             300x100
% Time:             0:0.2:500
% Relative speed    0 m/s

wavesStruct.seaState = 6;
wavesStruct.beta     = 3 * pi / 4;
wavesStruct.xVec     = linspace(0, 299, 300);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:200;
wavesStruct.U        = 0;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec, wavesStruct.U);
wavesStruct.waveType = 'long';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Wave #4: Sea state 3 wave coming from bow (front)
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       180 degrees (from the stern)
% Grid:             300x100
% Time:             0:0.2:500
% Relative speed    0 m/s

wavesStruct.seaState = 3;
wavesStruct.beta     = pi;
wavesStruct.xVec     = linspace(0, 299, 300);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:300;
wavesStruct.U        = 0;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec, wavesStruct.U);
wavesStruct.waveType = 'long';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Wave #5: Sea state 3 wave coming from bow (front), 15 knots
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): long-crested (unidirectional)
% Wave angle:       180 degrees (from the stern)
% Grid:             300x100
% Time:             0:0.2:500
% Relative speed    7.72 m/s

wavesStruct.seaState = 3;
wavesStruct.beta     = pi;
wavesStruct.xVec     = linspace(0, 299, 300);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:500;
wavesStruct.U        = 7.72;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec, wavesStruct.U);
wavesStruct.waveType = 'long';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Wave #6: Sea state 3 wave coming from bow (front), 15 knots
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): short-crested (multi-directional)
% Wave angle:       180 degrees (from the stern)
% Grid:             300x100
% Time:             0:0.2:500
% Relative speed    7.72 m/s

wavesStruct.seaState = 3;
wavesStruct.beta     = pi;
wavesStruct.xVec     = linspace(0, 299, 300);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:500;
wavesStruct.U        = 7.72;
dmu      = pi/20;
muVec    = -pi/2:dmu:pi/2;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec, ...
                                       wavesStruct.U,[],muVec,dmu);
wavesStruct.waveType = 'short';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Wave #7: Sea state 3 wave coming from bow (front)
% ------- Use wave with properties below
% Sea state:        3
% Wave type (beta): long-crested
% Wave angle:       180 degrees (from the stern)
% Grid:             2000x100
% Time:             0:0.2:50
% Relative speed    0 m/s

wavesStruct.seaState = 3;
wavesStruct.beta     = pi;
wavesStruct.xVec     = linspace(0, 849, 850);
wavesStruct.yVec     = linspace(0, 99, 100);
wavesStruct.Ts       = 0.2;
wavesStruct.tVec     = 0:wavesStruct.Ts:90;
wavesStruct.U        = 0;

wavesStruct.waves = simulateWaves(wavesStruct.seaState, ...
                                       wavesStruct.xVec, wavesStruct.yVec, ...
                                       wavesStruct.beta, wavesStruct.tVec, wavesStruct.U);
wavesStruct.waveType = 'long';

wavesStruct.displayName = [
   'Waves with properties: ', ...
   '\n   -Sea state: ', num2str(wavesStruct.seaState), ...
   '\n   -Significant wave height: ', num2str(getSignificantWaveHeight(wavesStruct.seaState)), ' m' ...
   '\n   -', wavesStruct.waveType, ' crested' ...
   '\n   -Main wave direction (beta): ', num2str(wavesStruct.beta), ' rad', ...
   '\n   -xVec: ', num2str(wavesStruct.xVec(1)), ':', num2str(1), ':', num2str(wavesStruct.xVec(end)), ' m', ...
   '\n   -yVec: ', num2str(wavesStruct.yVec(1)), ':', num2str(1), ':', num2str(wavesStruct.yVec(end)), ' m', ...
   '\n   -tVec: ', num2str(1), ':', num2str(wavesStruct.Ts), ':', num2str(wavesStruct.tVec(end)), ' s', ...
   '\n   -U: ', num2str(wavesStruct.U), ' m/s', ...
   '\n'];
saveWavesFile(wavesStruct);

%% Plot snapshot of short-crested seas for 3 different sea states
clear; clc;
beta     = 0;
xVec     = linspace(0, 300, 250);
yVec     = linspace(0, 300, 250);
Ts       = 0.2;
tVec     = 0:Ts:1;
dmu      = pi/2.5;
muVec    = -pi/2:dmu:pi/2;

waves_seaState_2 = simulateWaves(2, xVec, yVec, beta, tVec, 0, [], muVec, dmu);
waves_seaState_4 = simulateWaves(4, xVec, yVec, beta, tVec, 0, [], muVec, dmu);
% waves_seaState_6 = simulateWaves(5, xVec, yVec, beta, tVec, [], muVec, dmu);

figure(1);
subplot(1, 2, 1);
sea_1 = surf(xVec, yVec, waves_seaState_2(:, :, 1));
axis([10 250 10 250 -5 5]);
% title("Sea state 2");
xlabel("[m]"); ylabel("[m]"); zlabel("[m]");
subplot(1, 2, 2);
sea_3 = surf(xVec, yVec, waves_seaState_4(:, :, 1));
axis([10 250 10 250 -3.5 3.5]);
% title("Sea state 3");
xlabel("[m]"); ylabel("[m]"); zlabel("[m]");
% subplot(3, 1, 3);
% sea_6 = surf(xVec, yVec, waves_seaState_6(:, :, 1));
% title("Sea state 5");
% xlabel("[m]"); ylabel("[m]"); zlabel("[m]");
% axis([10 200 10 200 -4 4]);
                                   

%% Help file
function saveWavesFile(wavesStruct)
fileName = ['wave-files/waves__seaState_', num2str(wavesStruct.seaState), '__', ... 
            wavesStruct.waveType, ...
            '__beta_', num2str(round(wavesStruct.beta, 2)), ...
            '__grid_', num2str(length(wavesStruct.xVec)), 'x', num2str(length(wavesStruct.yVec)), ...
            '__time_0_', num2str(wavesStruct.Ts), '_', num2str(wavesStruct.tVec(end)), ...  
            '__U_', num2str(wavesStruct.U),...
            '.mat'];
save(fileName, 'wavesStruct');
disp(['Saved waves into file `', fileName, '`']);
end