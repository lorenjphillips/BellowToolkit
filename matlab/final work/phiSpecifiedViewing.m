function  phiSpecifiedViewing()

n = input('Enter the number of segments: ');
params = []; % Initialize

for i = 1:n
    fprintf('Enter data in [meters] for segment %d:\n', i);
    
    % Prompt for minimum and maximum values, and step size for the current segment
    min_l1 = input('Enter the MINimum value for l1: ');
    min_l2 = input('Enter the MINimum value for l2: ');
    min_l3 = input('Enter the MINimum value for l3: ');
    max_l1 = input('Enter the MAXimum value for l1: ');
    max_l2 = input('Enter the MAXimum value for l2: ');   
    max_l3 = input('Enter the MAXimum value for l3: ');
    d = input('Enter d: ');
    num_steps = input('Enter the number of steps: ');
    
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

% Prompt the user for the range of phi values they would like to analyze
min_phi_deg = input('Enter the MINimum value for phi (in degrees): ');
max_phi_deg = input('Enter the MAXimum value for phi (in degrees): '); % Could add validation

% Convert the range to radians
min_phi_rad = deg2rad(min_phi_deg);
max_phi_rad = deg2rad(max_phi_deg);

% Filter the params matrix based on the phi range
filtered_params = params(params(:,6) >= min_phi_rad & params(:,6) <= max_phi_rad, :);

% Display the filtered parameters
disp('Filtered Parameter Results (segment, l1, l2, l3, kappa, phi, ell):');
disp(filtered_params);

end

%{
g = robotindependentmapping(kappa, phi, ell, 20);
assignin('base','g',g)
    
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

% Display the result
disp(['The angle with respect to the x-y plane is ', num2str(angle_deg), ' degrees']);
disp(g)
end
%}
