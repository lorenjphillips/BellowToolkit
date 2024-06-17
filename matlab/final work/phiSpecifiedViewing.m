function  phiSpecifiedViewing()
% Prompt the user for the number of segments
n = input('Enter the number of segments: ');

% Initialize an empty matrix to store the results
results = [];

% Loop through each segment
for i = 1:n
    fprintf('Enter data in [meters] for segment %d:\n', i);
    
    % Prompt for minimum and maximum values of l1, l2, l3, and step size for the current segment
    min_l1 = input('Enter the minimum value for l1: ');
    max_l1 = input('Enter the maximum value for l1: ');
    min_l2 = input('Enter the minimum value for l2: ');
    max_l2 = input('Enter the maximum value for l2: ');
    min_l3 = input('Enter the minimum value for l3: ');
    max_l3 = input('Enter the maximum value for l3: ');
    num_steps = input('Enter the number of steps: ');
    d = input('Enter d: ');
    
       % Calculate the step increments based on the number of steps
    step_l1 = (max_l1 - min_l1) / (num_steps - 1);
    step_l2 = (max_l2 - min_l2) / (num_steps - 1);
    step_l3 = (max_l3 - min_l3) / (num_steps - 1);
    
    % Nested loops to iterate through all possible combinations of l1, l2, and l3 values
    for l1 = min_l1:step_l1:max_l1
        for l2 = min_l2:step_l2:max_l2
            for l3 = min_l3:step_l3:max_l3
                % Calculate kappa, phi, and ell
                kappa = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));
                phi = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));
                ell = (l1 + l2 + l3) / 3;
                
                % Store the segment number, current l1, l2, l3, kappa, phi, and ell in the results matrix
                results = [results; i, l1, l2, l3, kappa, phi, ell];
            end
        end
    end
end

% Display the results
disp('Results (segment, l1, l2, l3, kappa, phi, ell):');
disp(results);




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
