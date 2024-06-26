function  phiSpecifiedViewing()

n = input('Enter the number of segments: ');
params = []; % Initialize

for i = 1:n
    fprintf('Enter data in [meters] for segment %d:\n', i);
   
    min_l1 = .01; min_l2 = min_l1;  min_l3 = min_l1;  % Test Values
    max_l1 = .03; max_l2 = max_l1; max_l3 = max_l1 ;  % Test Values
    d = .005; num_steps = 7;  % Test Values
    
    %{
    % Prompt for minimum and maximum values, and step size for the current segment
        min_l1 = input('Enter the MINimum value for l: ');
        min_l2 = min_l1;%input('Enter the MINimum value for l2: ');
        min_l3 = min_l1;%input('Enter the MINimum value for l3: ');
        max_l1 = input('Enter the MAXimum value for l: ');
        max_l2 = max_l1;%input('Enter the MAXimum value for l2: ');   
        max_l3 = max_l1;%input('Enter the MAXimum value for l3: ');
        d = input('Enter d: ');
        num_steps = input('Enter the number of steps: ');
    %}

    step_l1 = (max_l1 - min_l1) / (num_steps - 1); % Calculate the step increments
    step_l2 = (max_l2 - min_l2) / (num_steps - 1);
    step_l3 = (max_l3 - min_l3) / (num_steps - 1);
 
    for l1 = min_l1:step_l1:max_l1   % Nested loops to iterate through all possible combinations.
        for l2 = min_l2:step_l2:max_l2
            for l3 = min_l3:step_l3:max_l3
                                     % Calculate kappa, phi, and ell
                kappa = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));
                phi = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));
                ell = (l1 + l2 + l3) / 3;
                
                params = [params; i, l1, l2, l3, kappa, phi, ell];
            end
        end
    end
end

disp('Check of Parameter Results (segment, l1, l2, l3, kappa, phi, ell):');
disp(size(params)); % Verify Dimensions %disp(params);

%% Running Robot Mapping Function

params_with_mapping = []; % Initialize a new matrix to store parameters with mapping indices
mapping_results = []; % Array to store all mapping results
angles_degrees = []; % Initialize an array to store the angles
for idx = 1:size(params, 1) % Iterate over 'params' and run each kappa, phi, and ell --> robotindependentmapping
    if params(idx, 1) == 1 % Only process if the segment is 1
        kappa = params(idx, 5);
        phi = params(idx, 6);
        ell = params(idx, 7);
        
        % Call the robotindependentmapping function and store the result
        mapping_result = robotindependentmapping(kappa, phi, ell, 10);
        
        % Store the mapping result in the array
        mapping_results = cat(3, mapping_results, mapping_result); % Build a 3D array. 
        % Append the parameters and the index to the new matrix
        params_with_mapping = [params_with_mapping; params(idx, :), idx];
        % Extract the components of the vector from the last row of mapping_result
        vx = mapping_result(end, 9);
        vy = mapping_result(end, 10);
        vz = mapping_result(end, 11);
        
        % Calculate the magnitude of the vector
        v_mag = sqrt(vx^2 + vy^2 + vz^2);
        
        % Calculate the magnitude of the projection of the vector onto the x-y plane
        v_proj_mag = sqrt(vx^2 + vy^2);
        
        % Check if the values are valid before computing the angle
        if v_mag > 0 && v_proj_mag / v_mag <= 1 && v_proj_mag / v_mag >= -1
            % Calculate the angle between the vector and the x-y plane
            angle_rad = acos(v_proj_mag / v_mag);
            
            % Convert the angle from radians to degrees
            angle_deg = rad2deg(angle_rad);
            
            % Append the angle to the array
            angles_degrees = [angles_degrees; angle_deg];
            
            % Append the angle to the params_with_mapping matrix
            params_with_mapping(end, 9) = angle_deg;
        else
            % If the values are not valid, append NaN
            angles_degrees = [angles_degrees; NaN];
            params_with_mapping(end, 9) = NaN;
        end
    end
end
    %disp('Updated Parameter Results (segment, l1, l2, l3, kappa, phi, ell, mapping_index):');
    %disp(params_with_mapping);

%% Plotting

% Prompt the user if they want to plot the results
plot_choice = input('Would you like to plot the results? (y/n): ', 's');

if strcmpi(plot_choice, 'y')
    % Plotting
    figure;
    fig=figure;
    fig.Color = [1 1 1];
    hold on;
    all_points = []; % Initialize variables for dynamic scaling
    for idx = 1:size(mapping_results, 3)
        g = mapping_results(:, :, idx);
        vx = g(end, 9); % Extract the components of the vector
        vy = g(end, 10);
        vz = g(end, 11);
        x = g(end, 13); % Extract the origin of the vector
        y = g(end, 14);
        z = g(end, 15);
        all_points = [all_points; x, y, z, vx, vy, vz];
    end
    scale_factor = 0.1 * max(max(all_points(:, 1:3)) - min(all_points(:, 1:3))); % Dynamic scale factor
    for i = 1:size(all_points, 1)
        quiver3(all_points(i, 1), all_points(i, 2), all_points(i, 3), ...
                all_points(i, 4), all_points(i, 5), all_points(i, 6), ...
                scale_factor, 'LineWidth', 3, 'Color', [0 0 1]);
    end

    hold off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Vector Plot');
    grid on;
end

%% Histogram
% Prompt the user if they want to plot the results
plot_choice2 = input('Would you like a histogram of the results? (y/n): ', 's');

if strcmpi(plot_choice2, 'y')
    % Display the updated params_with_mapping matrix
    disp('Updated Parameter Results (segment, l1, l2, l3, kappa, phi, ell, mapping_index, angle_degrees):');
    disp(params_with_mapping);
    
    % Plot the histogram of angles
    figure;
    histogram(angles_degrees, 'BinWidth', 5); % Adjust 'BinWidth' as needed
    xlabel('Angle (degrees)');
    ylabel('Frequency');
    title('Histogram of Angles');
    grid on;
end

end

%% Angle Calculations 1
%{
angles_degrees = []; % Calculate the angle for each mapping result and append to params_with_mapping
for idx = 1:size(mapping_results, 3)
    g = mapping_results(:, :, idx);
    vx = g(end, 9); % Extract the components of the vector
    vy = g(end, 10);
    vz = g(end, 11);
    v_mag = sqrt(vx^2 + vy^2 + vz^2); % Calculate the magnitude of the vector
    v_proj_mag = sqrt(vx^2 + vy^2); % Calculate the magnitude of the projection of the vector onto the x-y plane 
    angle_rad = acos(v_proj_mag / v_mag); % Calculate the angle between the vector and the x-y plane
    angle_deg = rad2deg(angle_rad); % Convert the angle from radians to degrees
    angles_degrees = [angles_degrees; angle_deg];
end
params_with_angles = [params_with_mapping, angles_degrees]; % Append the angles to the params_with_mapping matrix
disp('Updated Parameter Results (segment, l1, l2, l3, kappa, phi, ell, mapping_index, angle_degrees):');
% disp(params_with_angles); % Display the updated params with angles

% Create histogram of angles
figure;
histogram(angles_degrees, 'BinWidth', 5); % Adjust 'BinWidth' as needed
xlabel('Angle (degrees)');
ylabel('Frequency');
title('Histogram of Angles');
grid on;
%}


%% EXTRA NOTES
%{

% Calculate the magnitude of the vector
v_mag = sqrt(vx^2 + vy^2 + vz^2);

% Calculate the magnitude of the projection of the vector onto the x-y plane
v_proj_mag = sqrt(vx^2 + vy^2);

% Calculate the angle between the vector and the x-y plane
angle_rad = acos(v_proj_mag / v_mag);

% Convert the angle from radians to degrees
angle_deg = rad2deg(angle_rad);
%}
%% PLOT SETUP

%{
    fig=figure;
    fig.Color = [1 1 1]; 
    hold on
    set(fig,'Position',[0 0 1280 1024]);

    % Axes, Labels
    clearance = 0.03;
    axis([-(max(abs(g(:,13)))+clearance) max(abs(g(:,13)))+clearance ...
          -(max(abs(g(:,14)))+clearance) max(abs(g(:,14)))+clearance ...
          0 curvelength+clearance])
    ax=gca;
    ax.GridAlpha=0.3;

    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    grid on
    view([0.5 0.5 0.5])
    daspect([1 1 1])

    col = linspace(0.2,0.8,numseg);

%}

%% MORE NOTES
%{
    % THE FIRST 3 ARE POINT, LAST 3 ARE VECTORS, These specify the starting point (origin) of the arrow. 
quiver3(g(end,13),g(end,14),g(end,15),g(end,9),g(end,10),g(end,11),0.01,'LineWidth',3,'Color',[0 0 1]);

%}

% Now it is mapping everything but not stacking, need to choose phi
% tolerance and stack, leave line to make this an optional input
    % Any segment two within 5 degrees of segment 1 can be used as second/third
    % etc. let mapping function do stacking work





% Run g here with one segment and two segment the sort back into the
% matrix the vx/y/z and x/y/z

%{

% Prompt the user for the range of phi values they would like to analyze
min_phi_deg = input('Enter the MINimum value for phi (in degrees): ');
max_phi_deg = input('Enter the MAXimum value for phi (in degrees): '); % Could add validation

% Convert the range to radians
min_phi_rad = deg2rad(min_phi_deg);
max_phi_rad = deg2rad(max_phi_deg);

% Filter the params matrix based on the phi range
filtered_params = params(params(:,6) >= min_phi_rad & params(:,6) <= max_phi_rad, :);

% Display the filtered parameters
disp('Check Dimensions of Filtered Parameter Results (segment, l1, l2, l3, kappa, phi, ell):');
disp(size(filtered_params)); % Verify Dimensions

end
%}

% Calculate engle to horizatnal
% 
% PROMPT List all final angles in table
%
% PROMPT and box plot them (angles)

% PROMPT 3D plot (use quiver and position info from the matrix
