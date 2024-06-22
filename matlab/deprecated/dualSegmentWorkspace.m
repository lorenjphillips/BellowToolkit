function dualSegmentWorkspace(r, nb1, nb2, d1, d2, step, t)
    % Validate input parameters
    validateInputs(r, d1, d2, step);

    % Segment properties
    [l_1_min, l_1_max, l_2_min, l_2_max] = segmentLengths(r, nb1, nb2, t);

    % Initialize arrays for each segment
    segment1 = calculateSegment(l_1_min, l_1_max, d1, step);
    segment2_template = calculateSegment(l_2_min, l_2_max, d2, step);

    % Initialize variables to store combined results
    all_positions = [];
    all_phi = [];
    all_kappa = [];
    all_ell = [];

    % Iterate over each end position of segment 1
    for idx = 1:size(segment1.positions, 1)
        last_position = segment1.positions(idx, :);
        last_phi = segment1.phi(idx);
        last_kappa = segment1.kappa(idx);
        last_theta = last_kappa * segment1.ell(idx);

        % Rotation due to last segment's curvature and orientation
        R = rotationMatrix(last_phi, last_theta);

        % Create new instance of segment 2
        segment2 = segment2_template; % Copy template of segment2
        for i = 1:size(segment2.positions, 1)
            segment2.positions(i, :) = (R * segment2.positions(i, :).' + last_position.').';
        end

        % Combine segment 2 data to the total results
        all_positions = [all_positions; segment2.positions];
        all_phi = [all_phi; segment2.phi];
        all_kappa = [all_kappa; segment2.kappa];
        all_ell = [all_ell; segment2.ell];
    end

    % Plot results with the number of points in segment 1 for differentiation
    plotNewResults(all_positions, all_phi, all_kappa, all_ell, size(segment1.positions, 1));
end


function validateInputs(r, d1, d2, step)
    if d1 <= 0 || d2 <= 0 || r <= 0
        error('Diameters and radius must be positive values');
    end
    if step <= 0
        error('Step must be a positive value');
    end
      % Validate input parameters
    if d1 > d2
        disp('The bellows are larger at the base and decrease in size as they get towards the distal end.');
    end
    if d1 < d2
        disp('The bellows are smaller at the base and increase in size as they get towards the distal end.');
    end
    if d1 == d2
        disp('The bellows maintain the same diameter');
    end
    if step <= 0
        error('Step must be a positive value');
    end
    if d1 <= 0
        error('Diameter must be a positive value');
    end
    if r <= 0
        error('Diameter must be a positive value');
    end
end

function [l_1_min, l_1_max, l_2_min, l_2_max] = segmentLengths(r, nb1, nb2, t)
    % t = 1.5; % Thickness in mm
    l_1_min = t * nb1;
    l_1_max = nb1 * 2 * pi * r;
    l_2_min = t * nb2;
    l_2_max = nb2 * 2 * pi * r;
end

function segment = calculateSegment(l_min, l_max, d, step)
    l_range = l_min:step:l_max;
    num_elements = numel(l_range);
    phi = zeros(num_elements^3, 1);
    ell = zeros(num_elements^3, 1);
    kappa = zeros(num_elements^3, 1);
    positions = zeros(num_elements^3, 3);
    index = 1;
    for l1 = l_range
        for l2 = l_range
            for l3 = l_range
                [phi_val, ell_val, kappa_val] = continuumRobotKinematics(l1, l2, l3, d);
                %[phi_val, ell_val, kappa_val] = calculateKinematics(l1, l2, l3, d);
                phi(index) = phi_val;
                ell(index) = ell_val;
                kappa(index) = kappa_val;
                positions(index, :) = calculatePosition(ell_val, kappa_val, phi_val);
                index = index + 1;
            end
        end
    end
    segment = struct('positions', positions, 'phi', phi, 'kappa', kappa, 'ell', ell);
end

function R = rotationMatrix(phi, theta)
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R = Rz * Ry;
end

function [phi_val, ell_val, kappa_val] = continuumRobotKinematics(l1, l2, l3, d)

    % Calculate the orientation angle phi(q)
    phi_val = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));

    % Calculate the effective length ell(q)
    ell_val = (l1 + l2 + l3) / 3;

    % Calculate the curvature kappa(q)
    kappa_val = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));

end

function position = calculatePosition(ell, kappa, phi)
    if kappa ~= 0
        R = 1 / kappa; % Radius of curvature
        theta = kappa * ell; % Bend angle
        dist = R - (R * cos(theta)); % Distance from origin
        z = R * sin(theta); % along the bending direction
        x = dist * cos(phi); % lateral displacement
        y = dist * sin(phi); % depth displacement
    else
        x = ell * cos(phi); % If no curvature, position extends in straight line
        y = ell * sin(phi);
        z = 0; % No bending, hence Z remains constant
    end
    position = [x, y, z];
end


function plotNewResults(all_positions, all_phi, all_kappa, all_ell, num_points_segment1)
    figure(1);
    f = figure;
    f.Color = [1 1 1]; 
    clf;  % Clear the figure to ensure a fresh start
    hold on;

    % Plot positions for the first segment in black
    scatter3(all_positions(1:num_points_segment1,1), all_positions(1:num_points_segment1,2), all_positions(1:num_points_segment1,3), 36, 'k', 'filled');

    % Plot positions for the second segment with a colormap based on their z-values
    if size(all_positions, 1) > num_points_segment1
        scatter3(all_positions(num_points_segment1+1:end,1), all_positions(num_points_segment1+1:end,2), all_positions(num_points_segment1+1:end,3), 36, all_positions(num_points_segment1+1:end,3), 'filled');
        colormap(jet); % Apply colormap to the second segment scatter plot
        cbar = colorbar; % Create a colorbar
        cbar.Label.String = 'Height (Z) in mm'; % Set the label for the colorbar

        % Calculate direction vectors for each tip position of the second segment
        % Only add quivers at the endpoints of each instance of the second segment
        % Assuming each second segment has the same number of points
        num_points_per_instance = length(all_phi) / num_points_segment1;
        end_indices = num_points_segment1+1:num_points_per_instance:size(all_positions, 1);
        u = cos(all_phi(end_indices)) .* sin(all_kappa(end_indices) .* all_ell(end_indices));
        v = sin(all_phi(end_indices)) .* sin(all_kappa(end_indices) .* all_ell(end_indices));
        w = cos(all_kappa(end_indices) .* all_ell(end_indices));
        
        % Plotting vectors as arrows at the end points of each second segment
        quiver3(all_positions(end_indices,1), all_positions(end_indices,2), all_positions(end_indices,3), u, v, w, 'r', 'AutoScale', 'on', 'AutoScaleFactor', 2, 'LineWidth', 1.5);
    end

    % Labels and title
    title('Dual Segments: First Segment in Black, Second Segment Color Mapped');
    xlabel('X (mm) - Lateral', 'FontSize', 12);
    ylabel('Y (mm) - Longitudinal', 'FontSize', 12);
    zlabel('Z (mm) - Vertical', 'FontSize', 12);

    grid on;
    axis equal;
    hold off;
end

