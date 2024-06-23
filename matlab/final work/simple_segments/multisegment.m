function multisegment(max_theta)
    % Constants
    ell = 1.05;
    phi = 0;
    numSteps = max_theta + 1; % Number of steps for theta values
    theta_values = linspace(0, max_theta, numSteps); % Array of theta values
    ell_steps = 21; % Number of steps for ell1 and ell2 combinations (0:5:100 gives 21 steps)
    ell1_values = linspace(0, ell, ell_steps);
    ell2_values = ell - ell1_values;

    % Calculate total number of iterations for the 3D array
    totalIterations = numSteps * numSteps * ell_steps;

    % Display the total number of iterations and ask user if they wish to proceed
    fprintf('Total number of iterations: %d\n', totalIterations);
    proceed = input('Do you wish to proceed? (yes/no): ', 's');
    
    if strcmpi(proceed, 'no')
        disp('Operation cancelled by user.');
        return;
    end

    % Placeholder for output matrices
    output_3D_array = [];
    output_2D_matrix = zeros(totalIterations, 6); % Columns: index, ell1, ell2, kappa1, kappa2, theta

    index = 1; % Initialize index for storing results

    % Loop through ell1 and ell2 values
    for j = 1:ell_steps
        ell1 = ell1_values(j);
        ell2 = ell2_values(j);

        % Loop through theta1 and theta2 values
        for t1 = 1:numSteps
            theta1 = theta_values(t1);
            kappa1 = (theta1 * pi) / (180 * ell1); % Calculate kappa1 for theta1

            for t2 = 1:numSteps
                theta2 = theta_values(t2);
                kappa2 = (theta2 * pi) / (180 * ell2); % Calculate kappa2 for theta2

                % Call the function and get the result
                n_seg = 21; % Total number of points
                result = robotindependentmapping([kappa1 kappa2], [phi phi], [ell1 ell2], n_seg);

                % Append result to 3D array
                if isempty(output_3D_array)
                    [n, m] = size(result);
                    output_3D_array = zeros(n, m, totalIterations);
                end
                output_3D_array(:, :, index) = result;

                % Store the index, ell1, ell2, kappa1, kappa2, and theta
                output_2D_matrix(index, :) = [index, ell1, ell2, kappa1, kappa2, theta1 + theta2];

                % Increment index
                index = index + 1;
            end
        end
    end

    % Save or display the results as needed
    disp('Size of 3D Array of Results:');
    disp(size(output_3D_array));

    disp('2D Matrix of Indices, Ell1, Ell2, Kappa, Theta (Deg):');
    disp(size(output_2D_matrix));

%% Plotting
    % Constants for backbone plot
    col = lines(totalIterations); % Color array for the segments, using lines colormap
    seg_end = n_seg; % Number of points in each segment, as per your input to robotindependentmapping

    figure;
    fig = figure;
    fig.Color = [1 1 1];
    hold on;
    for idx = 1:size(output_3D_array, 3)
        g = output_3D_array(:, :, idx);
        vx = g(end, 9); % Extract the components of the vector
        vz = g(end, 11);
        x = g(end, 13); % Extract the origin of the vector
        z = g(end, 15);
        
        % Plot the vector using quiver
        quiver(x, z, vx, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        
        % Plot the backbone
        plot(g(1:seg_end, 13), g(1:seg_end, 15), 'LineWidth', 2, 'Color', col(idx, :)); % Project to XZ plane
    end

    % Set plot labels and title
    xlabel('X (arbitrary unit)');
    ylabel('Z (arbitrary unit)');
    title('2D Vectors and Backbones on the XZ Plane from the Origin');
    xlim([0, 1.1]); % Adjust these limits based on your data
    ylim([0, 1.1]); % Adjust these limits based on your data

    grid on;
    hold off;
    end
   
%% Deprecated
%{
    phi = 0;
    numSteps = max_theta + 1;
    kappa_max = ((max_theta * pi) / (180 * ell)); 
    kappa_values = linspace(0, kappa_max, numSteps);
    % disp(kappa_values);
    
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
        vy = g(end, 10);
        vz = g(end, 11);
        x = g(end, 13); % Extract the origin of the vector
        y = g(end, 14);
        z = g(end, 15);
        
        % Plot the vector using quiver3
        
        quiver(x, z, vx, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        % quiver3(x, y, z, vx, vy, vz, 'AutoScale', 'on', 'AutoScaleFactor', 0.05, 'MaxHeadSize', 0.01);
        % Plot the backbone
        plot(g(1:seg_end, 13), g(1:seg_end, 15), 'LineWidth', 2, 'Color', col(idx, :)); % Project to XZ plane

            % If there are multiple segments, you can use a similar loop -
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
    title('3D Vectors from output\_3D\_array');
    xlim([0, 1.1]); % Adjust these limits based on your data
    ylim([0, 1.1]); % Adjust these limits based on your data

    grid on;
    hold off;
end
%}