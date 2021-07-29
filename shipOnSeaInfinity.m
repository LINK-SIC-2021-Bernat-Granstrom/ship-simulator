%function shipOnSea(waves, tVec, Ts)
% SHIPONSEA Ship time simulation. The forces considered are only the
% hidrostatic force by the water. Ignores: coriolis effect, added mass,
% propels force. The variables in the motion model refer to the ship's
% center of mass. To calculate the inertia, the ship was assumed to be a
% solid cuboid. Assumes the initial angles phi, th, psi are zero, i.e. the
% local coordinate system is aligned with the global one.
%
% Inputs:
%   - waves: surface object containing waves
%   - Ts:    sampling time [s]
clear;
fprintf('\nStarted shipOnSea simulation!\n');
%disp('1) Loading waves...');
tic;
%waves = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').waves;
%beta = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').beta;
%xVec  = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').xVec;
%yVec  = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').yVec;
%tVec  = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').tVec;
%Ts    = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').Ts;
%displayName = load('waves_seaState_6_lambda_200_long_beta_1.57_grid_400x150.mat').displayName;
%disp('Waves loaded:');
%fprintf(displayName);

Ts = 0.5;
tVec = 0:Ts:5;

% Constants:
%   - Ax: cross-sectional area of ship along ship's x-axis [m^2]
%   - Ay: cross-sectional area of ship along ship's y-axis [m^2]
%   - Az: cross-sectional area of ship along ship's z-axis [m^2]
%   - M:  ship's mass [kg]
%   - Ix: ship's momentum of inertia along x axis
%   - Iy: ship's momentum of inertia along y axis
%   - Iz: ship's momentum of inertia along z axis
%   - kv: damping constant for ship's velocity
%   - kw: damping constant for ship's rotational velocity
%   - ro: water density [kg/m^3]
%   - g:  gravity 9.81 [m/s^2]
%   - l:  length of ship (along x-axis) [m]
%   - w:  width of ship (along y-axis) [m]
%   - h:  height of ship (along z-axis) [m]

% Non-tunable constants
len = 137; width = 15; height = 16;
ro = 997; g = 9.81;   

% Tunable constants
kv = 5; 
kw = 1e10; 
Ax = width * height; 
Ay = height * len; 
Az = len * width; 
M = 2.5e6; % 2 500 000 [kg]
I = M/12 * diag([width^2 + height^2, len^2 + height^2, len^2 + width^2]);
Ix = I(1, 1) * 1e-10; 
Iy = I(2, 2); 
Iz = I(3, 3);

% Create wave parameters
disp('1) Creating wave spectrum...');


beta = 0;
dw = 0.025;
w = (0.5:dw:3)';
wVec = zeros(length(w),1);
for k=1:length(w)        
    wVec(k) = w(k) - dw/2 + dw * rand;
end

dmu = pi/50;
mu = (-pi/2:dmu:pi/2);
muVec = zeros(1, length(mu));
for k=1:length(mu)        
    muVec(k) = mu(k) - dmu/2 + dmu * rand;
end

e = zeros(length(muVec), length(wVec));
for k=1:length(wVec) 
    for i=1:length(muVec)
        e(i, k) = 2 * pi * rand;
    end 
end

S = createSpectrum(3,wVec);
disp('Spectrum created!');

%seaState2waveHeightPoint(beta, t, muVec, dmu, wVec, dw, e, x, y, S)
% Read boat file (faces, vertices, normals)
disp('2) Loading stl boat file and computing normal vectors to faces...');
[F, V, N] = stlreadOwn('filteredNorfolkNew.stl');

% V = V + [300 300 4.25]; 
V = V + [30 30 4.25]; 
cog = [(max(V(:,1))+min(V(:,1)))/2,(max(V(:,2))+min(V(:,2)))/2,(max(V(:,3))+min(V(:,3)))/2];
cog(1) = cog(1) - 4.77;
cog(2) = cog(2) + 0.022;
cog(3) = cog(3) - 2;
j = 1;
for i=1:3:length(V)-2
    facePoint = [mean(V(i:i+2,1)), mean(V(i:i+2,2)), mean(V(i:i+2,3))];
    facePoints(j, :) = [mean(V(i:i+2,1)), mean(V(i:i+2,2)), mean(V(i:i+2,3))]; % The center of the face
    faceArea(j) = areaOfFace(V(i:i+2,1), V(i:i+2,2), V(i:i+2,3)); % The area of the face
    p0 = V(i, :)'; p1 = V(i+1,:)'; p2 = V(i+2,:)';
    n = cross(p0-p1, p0-p2)'./norm(cross(p0-p1, p0-p2)'); 
    normals(j, :) = -n; % The normal of the face
    j = j + 1;
end
[row, col] = find(isnan(normals));
facePoints(isnan(facePoints)) = 0;
normals(isnan(normals)) = 0;
disp('Boat loaded!');

% State: [x y z   u v w   phi th psi   w_phi w_th w_psi]'
%    x y z       u              v        w        phi th psi w_phi w_th w_psi
A = [0 0 0       1              0        0         0  0   0    0    0    0;
     0 0 0       0              1        0         0  0   0    0    0    0;
     0 0 0       0              0        1         0  0   0    0    0    0;
     0 0 0 -0.5*kv*ro*Ax/M      0        0         0  0   0    0    0    0;
     0 0 0       0   -0.5*kv*ro*Ay/M     0         0  0   0    0    0    0;
     0 0 0       0              0  -0.5*kv*ro*Az/M 0  0   0    0    0    0;
     0 0 0       0              0        0         0  0   0    1    0    0;
     0 0 0       0              0        0         0  0   0    0    1    0;
     0 0 0       0              0        0         0  0   0    0    0    1;
     0 0 0       0              0        0         0  0   0 -kw/(Ix*50)  0    0;
     0 0 0       0              0        0         0  0   0    0 -kw/Iy  0;
     0 0 0       0              0        0         0  0   0    0    0 -kw/Iz]; 
 
%    F_x/M   F_y/M   F_z/M   Tau_x/Ix   Tau_y/Iy   Tau_z/Iz
B = [  0       0       0        0          0          0;  % x dot
       0       0       0        0          0          0;  % y dot
       0       0       0        0          0          0;  % z dot
       1       0       0        0          0          0;  % u dot
       0       1       0        0          0          0;  % v dot
       0       0       1        0          0          0;  % w fot
       0       0       0        0          0          0;  % phi dot
       0       0       0        0          0          0;  % th dot
       0       0       0        0          0          0;  % psi dot
       0       0       0        1          0          0;  % w_phi dot
       0       0       0        0          1          0;  % w_th dot
       0       0       0        0          0          1]; % w_psi dot
C = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0];
sys = c2d(ss(A, B, C, []), Ts);

states = zeros(12, length(tVec));
% State:       [ xyz  u v w phi th psi   w_phi w_th w_psi]'
states(:, 1) = [ cog  0 0 0  0  0   0      0     0    0]';

% Rotate ship hull according to initial states
facePoints = (R(states(7,1), states(8,1), states(9,1))' * (facePoints - cog)')' + cog; 
normals = (R(states(7,1), states(8,1), states(9,1))' * normals')';
disp('3) Simulating states through time...');
% State update for all time steps
cogVec = zeros(length(tVec),3);
for tIdx=1:length(tVec)-1
    cogVec(tIdx, :) = cog;
    % ------- Compute sum of all forces F_net & sum of all torques Tau_net
    F_net = zeros(3, 1);
    Tau_net = zeros(3, 1);
    disp('before normals')
    for nIdx=1:length(normals)
        % Get index of waves that is closest to the evaluated normal
        xIdx = facePoints(nIdx, 1); 
        yIdx = facePoints(nIdx, 2); 
        
        facePointHeight = facePoints(nIdx, 3);
        eta = seaState2waveHeightPoint(beta, tIdx, muVec, dmu, wVec, dw, e, xIdx, yIdx, S);
        % Add nettoforce component if wave is above facePoint
        if eta > facePointHeight       
            % "Spring" force
            h = eta - facePointHeight;                                    
            volume = h * faceArea(nIdx);
            F_buoy = (ro * g * volume) * normals(nIdx, :)';
            F_net = F_net + F_buoy;      
                        
            % Torque
            lever = (facePoints(nIdx, :) - cog)'; % 3x1
            phi = states(7,tIdx); th = states(8,tIdx); psi = states(9,tIdx);
            Tau = cross(R(phi, th, psi) * lever, R(phi, th, psi) * F_buoy);
            Tau_net = Tau_net + Tau;
        end
    end
    disp('after normals')
    moment(:, tIdx) = Tau_net;
    forces(:, tIdx) = F_net;
    % ------- Set up input
    phi = states(7,tIdx); th = states(8,tIdx); psi = states(9,tIdx);
    forcesInLocalCoord = R(phi, th, psi) * [F_net(1)/M; F_net(2)/M; F_net(3)/M - g];
    torqueInLocalCoord = [Tau_net(1)/Ix; Tau_net(2)/Iy; Tau_net(3)/Iz];
    inputs = [forcesInLocalCoord' torqueInLocalCoord']';
    
    % ------- Time-update
    u = states(4,tIdx); v = states(5,tIdx); w = states(6,tIdx);
    w_phi = states(10,tIdx); w_th = states(11,tIdx); w_psi = states(12,tIdx);
    velInGlobalCoord = R(phi, th, psi)' * [u; v; w];
    rotVelDerivative = T(phi, th) * [w_phi; w_th; w_psi];
    % ------- Control (0 yaw):
%     inputs(6) = inputs(6)-1*states(9,tIdx);

    % Make time-update as below due to nonlinearities in the A matrix
    states(:, tIdx+1) = [states(1,tIdx) + sys.A(1,4) * velInGlobalCoord(1);
                         states(2,tIdx) + sys.A(2,5) * velInGlobalCoord(2);
                         states(3,tIdx) + sys.A(3,6) * velInGlobalCoord(3);
                         sys.A(4,4) * states(4,tIdx);
                         sys.A(5,5) * states(5,tIdx);
                         sys.A(6,6) * states(6,tIdx);
                         states(7,tIdx) + sys.A(7,10) * rotVelDerivative(1);
                         states(8,tIdx) + sys.A(8,11) * rotVelDerivative(2);
                         states(9,tIdx) + sys.A(9,12) * rotVelDerivative(3);
                         sys.A(10,10) * states(10,tIdx);
                         sys.A(11,11) * states(11,tIdx);
                         sys.A(12,12) * states(12,tIdx)] + sys.B * inputs;
    
    % ------- Update ship hull
    deltaX = states(1, tIdx+1) - states(1, tIdx);
    deltaY = states(2, tIdx+1) - states(2, tIdx);
    deltaZ = states(3, tIdx+1) - states(3, tIdx);
    
    deltaPhi = states(7, tIdx+1) - states(7, tIdx);
    deltaTh = states(8, tIdx+1) - states(8, tIdx);
    deltaPsi = states(9, tIdx+1) - states(9, tIdx);
    
    % Rotate facepoints and normals with R^t - correct?
    facePoints = (R(deltaPhi, deltaTh, deltaPsi)' * (facePoints - cog)')' + cog + [deltaX deltaY deltaZ]; 
    
    normals = (R(deltaPhi, deltaTh, deltaPsi)' * normals')';
    cog = cog + [deltaX deltaY deltaZ];
end
disp('Simulation done!');
toc;

%% Plot states
% Global position states x, y, z
figure;
subplot(4, 3, 1); plot(tVec * Ts, states(1, :)); title('Global CoG position X [m]'); xlabel('Time [s]');
subplot(4, 3, 2); plot(tVec * Ts, states(2, :)); title('Global CoG position Y [m]'); xlabel('Time [s]');
subplot(4, 3, 3); plot(tVec * Ts, states(3, :)); title('Global CoG position Z [m]'); xlabel('Time [s]');

% Local velocity states u, v, w
subplot(4, 3, 4); plot(tVec * Ts, states(4, :) / Ts); title('Local velocity U [m/s]'); xlabel('Time [s]');
subplot(4, 3, 5); plot(tVec * Ts, states(5, :) / Ts); title('Local velocity V [m/s]'); xlabel('Time [s]');
subplot(4, 3, 6); plot(tVec * Ts, states(6, :) / Ts); title('Local velocity W [m/s]'); xlabel('Time [s]');

% Local angles phi, th, psi
subplot(4, 3, 7); plot(tVec * Ts, rad2deg(states(7, :))); title('Roll \phi [deg]'); xlabel('Time [s]');
subplot(4, 3, 8); plot(tVec * Ts, rad2deg(states(8, :))); title('Pitch \theta [deg]'); xlabel('Time [s]');
subplot(4, 3, 9); plot(tVec * Ts, rad2deg(states(9, :))); title('Yaw \psi [deg]'); xlabel('Time [s]');

% Local rotational velocities w_phi, w_th, w_psi
subplot(4, 3, 10); plot(tVec * Ts, rad2deg(states(10, :)) / Ts); title('Roll velocity w_\phi [deg/s]'); xlabel('Time [s]');
subplot(4, 3, 11); plot(tVec * Ts, rad2deg(states(11, :)) / Ts); title('Pitch velocity w_\theta [deg/s]'); xlabel('Time [s]');
subplot(4, 3, 12); plot(tVec * Ts, rad2deg(states(12, :)) / Ts); title('Yaw velocity w_\psi [deg/s]'); xlabel('Time [s]');

% sgtitle(['Long crested from angle ', num2str(beta)]);

%% Simulate in 3D
figure;
p = patch('Faces',F,'Vertices',V,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.35); hold on;
axis('image');
view([-135 35]);
camlight('headlight');
material('dull');

% Sea with no waves:
% sea = surf(xVec, yVec, zeros(length(yVec), length(xVec)));
% Sea with waves:
sea = surf(xVec, yVec, waves(:, :, 1));
% axis([250 500 250 350 -20 20]);
axis([150 400 0 60 -20 20]);
xlabel('x [m]');
ylabel('y [m]');
%%
for t=2:length(tVec)-1
    sea.ZData = waves(:, :, t); % comment out for sea without waves
    deltaX = states(1, t) - states(1, t-1);
    deltaY = states(2, t) - states(2, t-1);
    deltaZ = states(3, t) - states(3, t-1);
    p.Vertices = (R(states(7,t), states(8,t), states(9,t))'*(V-cogVec(t,:))')'+[deltaX deltaY deltaZ] + cogVec(t,:);

    pause(0.001);
end

%% Help functions
function R = R(phi, th, psi)
% Rotation matrix using the xyz rule.
% Gustafsson, "Statistical Sensor Fusion" 3rd edition
% Eq. (13.7), p. 348
    R = [cos(th)*cos(psi)       cos(th)*sin(psi)       -sin(th);
         -cos(phi)*sin(psi)+sin(phi)*sin(th)*cos(psi) cos(phi)*cos(psi)+sin(phi)*sin(th)*sin(psi) sin(phi)*cos(th);
         sin(phi)*sin(psi)+cos(phi)*sin(th)*cos(psi) -sin(phi)*cos(psi)+cos(phi)*sin(th)*sin(psi) cos(phi)*cos(th)];
end
function T=T(phi, th)
% Gustafsson, "Statistical Sensor Fusion" 3rd edition
% Eq. (13.9), p. 349
    T = [1 sin(phi)*tan(th) cos(phi)*tan(th);
         0       cos(phi)         -sin(phi);
         0 sin(phi)/cos(th) cos(phi)/cos(th)];
     
    T = eye(3); % ok approximation as angle derivatives are small
end

%end