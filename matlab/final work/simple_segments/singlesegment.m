function singlesegment(max_theta)
    % Constants - EVENTUALLY INPUT ANGLE AND ELL
    ell = 1.05;
    phi = 0;
    numSteps = max_theta + 1;
    kappa_max = ((max_theta * pi) / (180 * ell)); % Check
    kappa_values = linspace(0, kappa_max, numSteps);
    disp(kappa_values);
    
    % Placeholder for output matrices
    output_3D_array = [];
    output_2D_matrix = zeros(numSteps, 3);
    
    % Loop through kappa values
    for i = 1:numSteps
        kappa = kappa_values(i);
        
        % Call the function and get the result
        n_seg=21;
        result = robotindependentmapping(kappa, phi, ell, n_seg);
        
        % Append result to 3D array
        if isempty(output_3D_array)
            [n, m] = size(result);
            output_3D_array = zeros(n, m, numSteps);
        end
        output_3D_array(:,:,i) = result;
        
        % Calculate viewing angle theta
        theta = (180 * kappa * ell) / pi;
        
        % Store the index, kappa value, and theta value
        output_2D_matrix(i, :) = [i, kappa, theta];
    end
    
    % Save or display the results as needed
    disp('3D Array of Results:');
    % disp(output_3D_array);
    
    disp('2D Matrix of Indices, Kappa Values, Theta (Deg):');
    disp(output_2D_matrix);

%% Plotting vectors
   % Constants for backbone plot
col = lines(numSteps); % Color array for the segments, using lines colormap
seg_end = n_seg; % Number of points in each segment, as per your input to robotindependentmapping


figure;
    fig = figure;
    fig.Color = [1 1 1];
    hold on;
    for idx = 1:size(output_3D_array, 3)
        g = output_3D_array(:, :, idx);
        vx = g(end, 9); % Extract the components of the vector
        % vy = g(end, 10);
        vz = g(end, 11);
        x = g(end, 13); % Extract the origin of the vector
        % y = g(end, 14);
        z = g(end, 15);
        
        % Plot the vector using quiver3
        
        quiver(x, z, vx, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        % quiver3(x, y, z, vx, vy, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        % Plot the backbone
        plot(g(1:seg_end, 13), g(1:seg_end, 15), 'LineWidth', 2, 'Color', col(idx, :)); % Project to XZ plane

            %  Formultiple segments,  use a similar loop -
            % will need to add numseg abovve (initilize and answer)
    % for i = 2:numseg
    %     plot(g(seg_end(i-1):seg_end(i), 13), g(seg_end(i-1):seg_end(i), 15), ...
    %          'LineWidth', 2, 'Color', col(idx, :));
    % end
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