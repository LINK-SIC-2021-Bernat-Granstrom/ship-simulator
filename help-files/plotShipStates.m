function plotShipStates(states, tVec, Ts, beta)
figure;
    % Ignore first 10 time steps as the system still is stabilizing itself
    tVec = tVec(10:end);
    states = states(:, 10:end);
    
    % Global position states x, y, z
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

    sgtitle(['Long crested from ', num2str(rad2deg(beta)), ' degrees']);
end