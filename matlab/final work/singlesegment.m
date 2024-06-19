function singlesegment()
    % Constants - EVENTUALLY INPUT ANGLE AND ELL
    ell = 1;
    phi = 0;
    numSteps = 91;
    kappa_max = (9 * ell) / (5 * pi);
    kappa_values = linspace(0, kappa_max, numSteps);
    disp(kappa_values);
    
    % Placeholder for output matrices
    output_3D_array = [];
    output_2D_matrix = zeros(numSteps, 2);
    
    % Loop through kappa values
    for i = 1:numSteps
        kappa = kappa_values(i);
        
        % Call the function and get the result
        result = robotindependentmapping(kappa, phi, ell, 20);
        
        % Append result to 3D array
        if isempty(output_3D_array)
            [n, m] = size(result);
            output_3D_array = zeros(n, m, numSteps);
        end
        output_3D_array(:,:,i) = result;
        
        % Store the index and kappa value
        output_2D_matrix(i, :) = [i, kappa];
    end
    
    % Save or display the results as needed
    disp('3D Array of Results:');
    % disp(output_3D_array);
    
    disp('2D Matrix of Indices and Kappa Values:');
    disp(output_2D_matrix);


% Plotting
figure;
fig=figure;
fig.Color = [1 1 1];
hold on;
all_points = []; % Initialize variables for dynamic scaling
for idx = 1:size(output_3D_array, 3)
    g = output_3D_array(:, :, idx);
    vx = g(end, 9); % Extract the components of the vector
    vy = g(end, 10);
    vz = g(end, 11);
    x = g(end, 13); % Extract the origin of the vector
    y = g(end, 14);
    z = g(end, 15);
    all_points = [all_points; x, y, z, vx, vy, vz];
end
for i = 1:size(all_points, 1)
    quiver3(all_points(i, 1), all_points(i, 2), all_points(i, 3), ...
            all_points(i, 4), all_points(i, 5), all_points(i, 6), ...
            0.05, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', .02);

end

hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Plot');
grid on;
end
