# Ship Simulator
Ship simulator is made of MATLAB files that simulates the **motion of waves** and the **motion of ships** on those waves. It was developed as part of a 6-week summer project between the competence center [LINK-SIC](https://liu.se/en/research/link-sic) and a company. 

For a quick overview of the results, see videos of a [pitch test](https://www.youtube.com/watch?v=bgkPYYFX2YQ) and a [roll test](https://www.youtube.com/watch?v=zRmFIrhpuu0).

## Installation
Just clone the repository on your machine. Only official MATLAB toolboxes are used.

## Demos
### Simulate waves
The demo-file for wave simulation is _demoCreateWaves.m_. See below how to simulate a wave of sea state 6.

```
%% Wave  #1:  Sea  state 6 wave  coming  from  bow (front)
% ------- Use  wave  with  properties  below
% Sea  state:          6
% Wave  type (beta): long-crested (unidirectional)
% Wave  angle:         180  degrees
% Grid:                200 x100
% Time:                0:0.2:500
% Relative  speed     0 m/s
wavesStruct.seaState = 6;
wavesStruct.beta      = pi;
wavesStruct.xVec      = linspace(0, 299, 300);
wavesStruct.yVec      = linspace(0, 99, 100);
wavesStruct.Ts         = 0.2;
wavesStruct.tVec      = 0: wavesStruct.Ts :500;
wavesStruct.U          = 0;
wavesStruct.waves = simulateWaves(wavesStruct.seaState, wavesStruct.xVec, wavesStruct.yVec, wavesStruct.beta, wavesStruct.tVec, wavesStruct.U);
wavesStruct.waveType ='long';
wavesStruct.displayName = ['Waves  with  properties:', '\n    -Sea  state:', num2str(wavesStruct.seaState), ...
                           '\n    -Significant  wave  height:', num2str(getSignificantWaveHeight(wavesStruct.seaState)),'m'...
                           '\n    -', wavesStruct.waveType ,'crested'...
                           '\n    -Main  wave  direction (beta):', num2str(wavesStruct.beta),'rad', ...
                           '\n    -xVec:', num2str(wavesStruct.xVec (1)),':', num2str (1),':', num2str(wavesStruct.xVec(end)),'m', ...
                           '\n    -yVec:', num2str(wavesStruct.yVec (1)),':', num2str (1),':', num2str(wavesStruct.yVec(end)),'m', ...
                           '\n    -tVec:', num2str (1),':', num2str(wavesStruct.Ts),':', num2str(wavesStruct.tVec(end)),'s', ...
                           '\n    -U:', num2str(wavesStruct.U),'m/s', ...
                           '\n'];
saveWavesFile(wavesStruct);
```

See an example of multi and unidirectional seas of sea state 2 and 4.

![Image of short and long crested](https://github.com/LINK-SIC-2021-Bernat-Granstrom/ship-simulator/blob/main/img/short-long-poster-cropped.png)

### Simulate ship on waves
The demo file for simulating a ship is _demoSimulateShip.m_. See below an example of how to simulate a ship on waves of sea state 6 created above.

```
%% ------- Create  ship  with  properties  below
% Ship:                            HMS  Norfolk (hull  only)
% Mass:                            2.5e
% Dimensions:                     137 x15x16
% Vertices  initial  position:   [30 35  5.55]
% Reference  velocity  along u: 0 [m/s] (should  always  be 0, the  ship  speed  is set in the  wave  file)
% Reference  yaw:                 0 degrees
% CoG  offset:                     [ -4.77  0.022  -2]
% Position  of  heliPad           [40 35 6]
shipStruct.file          ='filteredNorfolkNew.stl';
shipStruct.M             = 2.5e6;13shipStruct.len           = 137;
shipStruct.width         = 15;15shipStruct.height       = 16;
shipStruct.verticesPos = [30 35  5.55];17shipStruct.refSpeedU    = 0;
shipStruct.refYaw       = 0;19shipStruct.cogOffset    = [ -4.77  0.022  -2];
shipStruct.helipadPos   = [40 35 6];
%                  [v_u  v_v  v_w  phi th psi  w_phi  w_th  w_psi]'
shipStruct.x0 = [0     0    0   0    0   0     0     0     0]; 
% Initial  state  values2324%% --------------- Demo  #1:  Sea  state 6 pitch  test (waves  from  the  north)
% ------- Use  wave  with  properties  below
% Sea  state:          6
% Wave  type (beta): long -crested (unidirectional)
% Wave  angle:         180  degrees (from  the  front)
% Grid:                300 x100
% Time:                0:0.2:500
% Ship  speed          0
isPlot = true;
isVisual = true;
waveFile ='waves_seaState_6_long_beta_3 .14 _grid_300x100_time_0_0 .2 _500_U_0.mat';
[states , face , vert , cogVec] = simulateShip(waveFile , shipStruct , isPlot , isVisual);
```

A snapshot of the resulting simulation is shown below.

![Ship simulation](https://github.com/LINK-SIC-2021-Bernat-Granstrom/ship-simulator/blob/main/img/ship-on-wave-sea-state-6.png)


## License
[MIT](https://github.com/LINK-SIC-2021-Bernat-Granstrom/ship-simulator/blob/main/LICENSE)
