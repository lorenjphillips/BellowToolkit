function singlesegment(max_theta)
    % Constants - EVENTUALLY INPUT ANGLE AND ELL
    ell = 1.05; % Manually selecting effective length (arc length)
    phi = 0; % Manually selecting orientation of 0 degrees
    numSteps = max_theta + 1;
    kappa_max = ((max_theta * pi) / (180 * ell)); % Convert max theta to kappa
    kappa_values = linspace(0, kappa_max, numSteps);
    disp(kappa_values);
    
    % Initializing output matrices
    output_3D_array = [];
    output_2D_matrix = zeros(numSteps, 3);
    
    % Loop through kappa (curvature) values
    for i = 1:numSteps
        kappa = kappa_values(i);
        
        % Call robot mapping to get transformation matricies
        n_points=20; % Manually choosing 20 points per segment 
        result = robotindependentmapping(kappa, phi, ell, n_points);
        
        % Append result to 3D array
        if isempty(output_3D_array)
            [n, m] = size(result);
            output_3D_array = zeros(n, m, numSteps);
        end
        output_3D_array(:,:,i) = result;    
        theta = (180 * kappa * ell) / pi; % Calculate viewing angle theta
        output_2D_matrix(i, :) = [i, kappa, theta]; % Store the index, kappa value, and theta value
    end
    % disp('3D Array of Results:');
    % disp(output_3D_array);
    
    disp('2D Matrix of Indices, Kappa Values, Theta (Deg):');
    disp(output_2D_matrix);

%% Plotting vectors
   
col = lines(numSteps); % Made a color array for the segments, using lines colormap
seg_end = n_points; % Number of points in each segment, as per input to robotindependentmapping

figure;
    fig = figure;
    fig.Color = [1 1 1]; % I like to set the background to white
    hold on;
    for idx = 1:size(output_3D_array, 3)
        g = output_3D_array(:, :, idx);
        vx = g(end, 9); % Extract the components of the vector
        vy = g(end, 10);
        vz = g(end, 11);
        x = g(end, 13); % Extract the origin of the vector
        y = g(end, 14);
        z = g(end, 15);
        
        % Plot the vector using quiver3
        quiver(x, z, vx, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        % FOR 3D CASE...
        % quiver3(x, y, z, vx, vy, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01); 
        
        % Plot the backbone
        plot(g(1:seg_end, 13), g(1:seg_end, 15), 'LineWidth', 2, 'Color', col(idx, :)); % Project to XZ plane

            % When there are multiple segments, use a similar loop -
            % will need to add numseg abovve (initilize and answer)
                   %for i = 2:numseg
                     %  plot(g(seg_end(i-1):seg_end(i), 13), g(seg_end(i-1):seg_end(i), 15), ...
                     % 'LineWidth', 2, 'Color', col(idx, :));

end
    
    % Set plot labels and title
    xlabel('X (arbitrary unit)');
    ylabel('Z (arbitrary unit)');
    zlabel('Z');
    title('3D Vectors from output\_3D\_array');
    xlim([0, 1.1]); % Adjust these limits based on ell
    ylim([0, 1.1]); % Adjust these limits based on ell
        % Eventually - set above limits to be auto and keep scale
        % constant

    grid on;
    hold off;
end