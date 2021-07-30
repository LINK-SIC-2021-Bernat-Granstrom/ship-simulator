function visualizeSimulation(states, waves, xVec, yVec, tVec, faces, vertices, cogVec)
% VISUALIZESIMULATION Visualize 3D simulation of a ship states on given 
% waves.
    figure;
    p = patch('Faces', faces,'Vertices', vertices, 'FaceColor', [0.8 0.8 1.0], ...
             'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
             'AmbientStrength', 0.35); hold on;
    axis('image');
    view([-135 35]);
    camlight('headlight');
    material('dull');

    % ---------- Sea with no waves:
    % sea = surf(xVec, yVec, zeros(length(yVec), length(xVec)));
    % ---------- Sea with waves:
    sea = surf(xVec, yVec, waves(:, :, 1));
    
    % axis([250 500 250 350 -20 20]);
    axis([5 200 0 70 -20 20]);
    xlabel('x [m]');
    ylabel('y [m]');
    pause(0.1);
    deltaX = 0; deltaY = 0; deltaZ = 0;

    for t=2:length(tVec)-1
        sea.ZData = waves(:, :, t); % comment out for sea without waves
       
        deltaX = deltaX + states(1, t) - states(1, t-1);
        deltaY = deltaY + states(2, t) - states(2, t-1);
        deltaZ = deltaZ + states(3, t) - states(3, t-1);

        p.Vertices = (R(states(7,t), states(8,t), states(9,t))' * (vertices-cogVec(t,:))')'...
                     + cogVec(t,:) + [deltaX deltaY deltaZ];
        drawnow;
    end
end
