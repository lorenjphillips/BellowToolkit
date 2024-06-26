function ss = singlesegment(max_theta)
    ell = 1.00; % Constants - eventually to become inputs
    phi = 0;
    numSteps = max_theta + 1;
    kappa_max = ((max_theta * pi) / (180 * ell)); % Check
    kappa_values = linspace(0, kappa_max, numSteps);
    disp('Check of Kappa: '); disp(kappa_values);
    output_3D_array = []; % Init output matrices
    output_2D_matrix = zeros(numSteps, 5); % Init output matrices
    
    % Loop through kappa values
    for i = 1:numSteps
        kappa = kappa_values(i);
        n_seg=20;
        result = robotindependentmapping(kappa, phi, ell, n_seg); % Call the function
        
         % Append result to 3D array
        if isempty(output_3D_array)
            [n, m] = size(result);
            output_3D_array = zeros(n, m, numSteps);
        end
        output_3D_array(:,:,i) = result;
        theta = (180 * kappa * ell) / pi; % Calculate viewing angle theta

        g = result;  % Extract the components of the vector from the last row of g
        vx = g(end, 9); vy = g(end, 10); vz = g(end, 11);
        v_mag = sqrt(vx^2 + vy^2 + vz^2); % Calculate the magnitude of the vector
        v_proj_mag = sqrt(vx^2 + vy^2); % Calculate the magnitude of the projection of the vector onto the x-y plane
        angle_rad = acos(v_proj_mag / v_mag); % Calculate the angle between the vector and the x-y plane
        angle_deg = rad2deg(angle_rad); % Convert the angle from radians to degrees
        distance = g(end, 13); % Extract the origin of the vector
        output_2D_matrix(i, :) = [i, kappa, theta, angle_deg, distance]; % Store the index, kappa value, theta value, and calculated angle
    end
    
    % disp('3D Array of Results:'); % Display the results as needed
    % disp(output_3D_array);
    % disp('2D Matrix of Indices, Kappa Values, Theta Input (Deg), Theta Output wrt x-axis:');
     disp(output_2D_matrix);

     ss = output_2D_matrix(:, [4, 5]);

%% Plotting vectors
col = lines(numSteps); % Color array for the segments, using lines colormap
seg_end = n_seg; % Number of points in each segment, as per your input to robotindependentmapping

figure;
    fig = figure; fig.Color = [1 1 1]; hold on;
    for idx = 1:size(output_3D_array, 3)
        g = output_3D_array(:, :, idx);
        vx = g(end, 9); % Extract components of the vector
        vz = g(end, 11);
        x = g(end, 13); % Extract origin of the vector
        z = g(end, 15);
        
        quiver(x, z, vx, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01); % Plot the vector using quiver3
        % quiver3(x, y, z, vx, vy, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        plot(g(1:seg_end, 13), g(1:seg_end, 15), 'LineWidth', 2, 'Color', col(idx, :)); %  % Plot the backbone, Project to XZ plane

end
    
    % Set plot labels and title
    xlabel('X (arbitrary unit)');
    ylabel('Z (arbitrary unit)');
    zlabel('Z');
    title('Curve Backbones and Directions');
    xlim([0, 1.1]); 
    ylim([0, 1.1]); 

    grid on;
    hold off;
end