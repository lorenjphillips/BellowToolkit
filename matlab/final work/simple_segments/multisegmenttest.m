function mst = multisegmentloop(max_theta,ell)
    
    ell = 1.00; % Define total length of both segments
    phi = 0; % Orientation set to 0 for 2D curves
    numSteps = max_theta + 1; % Number of steps for theta values
    theta_values = linspace(0, max_theta, numSteps); % Array of theta values
    ell_steps = 21; % Number of steps for ell1 and ell2 combinations (0:5:100 gives 21 steps)
    ell1_values = linspace(0, ell, ell_steps); % Create array of possible segment 1 lengths
    ell2_values = ell - ell1_values; % Create array of possible segment 2 lengths
  
    totalIterations = numSteps * numSteps * ell_steps; % Calculate total number of iterations for the 3D array
    % Display arc lengths
    disp('ell1 and ell2 values:');
    disp(table(ell1_values', ell2_values', 'VariableNames', {'Seg 1 Lengths', 'Seg 2 Lengths'}));
        
    fprintf(2,'Total number of iterations: %d\n', totalIterations);
    proceed = input('Do you wish to proceed? (yes/no): ', 's');
    
    if strcmpi(proceed, 'no')
        disp('Operation cancelled by user.');
        return;
    end
    output_3D_array = []; % Placeholder for output matrices
    output_2D_matrix = zeros(totalIterations, 7); % Columns: index, ell1, ell2, kappa1, kappa2, theta1, theta2

    index = 1; % Initialize index for storing results

    for j = 1:ell_steps % Loop through ell1 and ell2 values (arc lengths)
        ell1 = ell1_values(j);
        ell2 = ell2_values(j);

        for t1 = 1:numSteps % Loop through theta1 and theta2 values - possible curvatures.
            theta1 = theta_values(t1);
            kappa1 = (theta1 * pi) / (180 * ell1); % Calculate kappa1 for theta1

            for t2 = 1:numSteps
                theta2 = theta_values(t2);
                kappa2 = (theta2 * pi) / (180 * ell2); % Calculate kappa2 for theta2

                % Call the mapping function from CRVis and save the result
                n_seg = 20; % Total number of points
                result = robotindependentmapping([kappa1 kappa2], [phi phi], [ell1 ell2], n_seg);

                if isempty(output_3D_array) % Append result to 3D array
                    [n, m] = size(result);
                    output_3D_array = zeros(n, m, totalIterations);
                end
                output_3D_array(:, :, index) = result;

                % Store the index, ell1, ell2, kappa1, kappa2, and theta
                output_2D_matrix(index, :) = [index, ell1, ell2, kappa1, kappa2, theta1, theta2];
                

                index = index + 1; % Increment index
            end
        end
    end

    % Prompt the user for the ell ratio
    selected_ratio = input('Enter the desired ell ratio as a percentage (e.g., 5 for 5:95): ');
    
    % Validate input
    if selected_ratio < 0 || selected_ratio > 100
        error('Invalid ratio. Please enter a value between 0 and 100.');
    end

    % Calculate the corresponding ell1 and ell2 values
    ell1_selected = (selected_ratio / 100) * ell;
    ell2_selected = ell - ell1_selected;

    % Find the indices corresponding to the selected ratio
    selected_indices = find(output_2D_matrix(:, 2) == ell1_selected & output_2D_matrix(:, 3) == ell2_selected);
    
    % Create the (nx2) matrix
    selected_matrix = zeros(length(selected_indices), 2);
    for idx = 1:length(selected_indices)
        index = selected_indices(idx);
        g = output_3D_array(:, :, index);
        total_theta = output_2D_matrix(index, 6) + output_2D_matrix(index, 7);
        x = g(end-1, 13);
        selected_matrix(idx, :) = [total_theta, x];
    end

    % Display the selected matrix
    disp('Selected Matrix (Total Theta, X Component):');
    fprintf('%.3f %.3f\n', selected_matrix.');

    %% Plotting
    % Constants for backbone plot
    col = lines(length(selected_indices)); % Color array for the segments, using lines colormap
    seg_end = (n_seg)*2; % Number of points in each segment, as per your input to robotindependentmapping

    figure;
    fig = figure;
    fig.Color = [1 1 1];
    hold on;
    for idx = 1:length(selected_indices)
        index = selected_indices(idx);
        g = output_3D_array(:, :, index);
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
    title(['2D Vectors and Backbones on the XZ Plane for ell1:ell2 = ' num2str(selected_ratio) ':' num2str(100 - selected_ratio)]);
    xlim([0, 1.1]); % Adjust these limits based on your data
    ylim([0, 1.1]); % Adjust these limits based on your data

    grid on;
    hold off;

    % Plot histogram for final theta values
    figure;
    histogram(selected_matrix(:, 1), 'BinWidth', 5); % Adjust BinWidth as needed
    xlabel('Total Theta (degrees)');
    ylabel('Frequency');
    title(['Histogram of Total Theta Values for ell1:ell2 = ' num2str(selected_ratio) ':' num2str(100 - selected_ratio)]);
    
    % Plot histogram for x components
    figure;
    histogram(selected_matrix(:, 2), 'BinWidth', 0.05); % Adjust BinWidth as needed
    xlabel('X Component (arbitrary unit)');
    ylabel('Frequency');
    title(['Histogram of X Components for ell1:ell2 = ' num2str(selected_ratio) ':' num2str(100 - selected_ratio)]);

    mst = selected_matrix; % Return the selected matrix as the function output
end
