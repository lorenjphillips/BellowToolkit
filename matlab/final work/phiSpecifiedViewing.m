function  phiSpecifiedViewing()

n = input('Enter the number of segments: ');
params = []; % Initialize

for i = 1:n
    fprintf('Enter data in [meters] for segment %d:\n', i);
    
    % Test Values
    min_l1 = 1; min_l2 = min_l1;  min_l3 = min_l1;
    max_l1 = 3; max_l2 = max_l1; max_l3 = max_l1 ;
    d = 5; num_steps = 3;
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
    % Calculate the step increments
    step_l1 = (max_l1 - min_l1) / (num_steps - 1);
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
disp(size(params)); % Verify Dimensions
disp(params);

% Initialize a new matrix to store parameters with mapping indices
params_with_mapping = [];

% Iterate over the params matrix and run each set of kappa, phi, and ell through robotindependentmapping
for idx = 1:size(params, 1)
        if params(idx, 1) == 1 
        kappa = params(idx, 5);
        phi = params(idx, 6);
        ell = params(idx, 7);
        
        % Call the robotindependentmapping function and store the result
        mapping_result = robotindependentmapping(kappa, phi, ell, 10);
        
        % Create a variable name dynamically to store each mapping result
        mapping_name = ['mapping_' num2str(idx)];
        eval([mapping_name ' = mapping_result;']);
        
        % Append the parameters and the index to the new matrix
        params_with_mapping = [params_with_mapping; params(idx, :), idx];
        
        end
end

% Display the updated params with mapping indices
disp('Updated Parameter Results (segment, l1, l2, l3, kappa, phi, ell, mapping_index):');
disp(params_with_mapping);
disp(size(params_with_mapping))



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

%{
g = robotindependentmapping(kappa, phi, ell, 20);
    
% Extract the components of the vector from the last row of g
vx = g(end, 9);
vy = g(end, 10);
vz = g(end, 11);

% Calculate the magnitude of the vector
v_mag = sqrt(vx^2 + vy^2 + vz^2);

% Calculate the magnitude of the projection of the vector onto the x-y plane
v_proj_mag = sqrt(vx^2 + vy^2);

% Calculate the angle between the vector and the x-y plane
angle_rad = acos(v_proj_mag / v_mag);

% Convert the angle from radians to degrees
angle_deg = rad2deg(angle_rad);

end
%}
