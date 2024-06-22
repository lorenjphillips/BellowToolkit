function plotFilteredContinuumVectors(l_min, l_max, d, step, phi_min, phi_max)
    % Validate input parameters
    validateInputs(l_min, l_max, d, step);

    % Initialize the range for each tendon
    l_range = l_min:step:l_max;
    num_elements = numel(l_range);
    
    % Initialize arrays
    phi = zeros(num_elements^3, 1);
    ell = zeros(num_elements^3, 1);
    kappa = zeros(num_elements^3, 1);
    positions = zeros(num_elements^3, 3);
    
    % Iterate over each tendon length
    index = 1;
    for l1 = l_range
        for l2 = l_range
            for l3 = l_range
                % Kinematics calculations
                [current_phi, current_ell, current_kappa] = calculateKinematics(l1, l2, l3, d);
                phi(index) = current_phi;
                ell(index) = current_ell;
                kappa(index) = current_kappa;
                positions(index, :) = calculatePosition(current_ell, current_kappa, current_phi);
                index = index + 1;
            end
        end
    end

    % Filter results based on phi range and plot
    indices = (phi >= phi_min) & (phi <= phi_max);
    if any(indices)
        filtered_positions = positions(indices, :);
        filtered_phi = phi(indices);
        filtered_kappa = kappa(indices);
        filtered_ell = ell(indices);

        plotResults(filtered_positions, filtered_phi, filtered_kappa, filtered_ell, l_min, l_max, d, step, phi_min, phi_max);
    else
        disp('No data points within the specified phi range.');
    end
end

function validateInputs(l_min, l_max, d, step)
    if l_min >= l_max
        error('Minimum tendon length must be less than maximum tendon length');
    end
    if step <= 0
        error('Step must be a positive value');
    end
    if d <= 0
        error('Diameter must be a positive value');
    end
end

function [phi, ell, kappa] = calculateKinematics(l1, l2, l3, d)
    phi = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));
    ell = (l1 + l2 + l3) / 3;
    kappa = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));
end

function position = calculatePosition(ell, kappa, phi)
    if kappa ~= 0
        R = 1 / kappa;
        theta = kappa * ell;
        
        dist = R - (R * cos(theta)); % Distance from origin
        z = R * sin(theta); % along the bending direction
        x = dist * cos(phi); % lateral displacement
        y = dist * sin(phi); % depth displacement
    else
        x = 0;
        y = 0;
        z = ell;
    end
    position = [x, y, z];
end

function plotResults(positions, phi, kappa, ell, l_min, l_max, d, step, phi_min, phi_max)
    figure;
    f = figure;
    f.Color = [1 1 1];              % Set the color to white
    [X, Y, Z] = cylinder(0.2, 50);
    surf(X, Y, Z, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    hold on;
    scatter3(positions(:,1), positions(:,2), positions(:,3), 'filled');

    % Calculate direction vectors
    num_points = size(positions, 1);
    u = zeros(num_points, 1);
    v = zeros(num_points, 1);
    w = zeros(num_points, 1);
    for i = 1:num_points
        theta = kappa(i) * ell(i);
        u(i) = cos(phi(i)) * sin(theta);
        v(i) = sin(phi(i)) * sin(theta);
        w(i) = cos(theta);
    end
    quiver3(positions(:,1), positions(:,2), positions(:,3), u, v, w, 0.5, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 1.5);

    title('Filtered Tip Positions and Directions of the Continuum Robot');
    xlabel('X (mm) - Lateral');
    ylabel('Y (mm) - Longitudinal');
    zlabel('Z (mm) - Vertical');
    grid on;
    axis equal;
    hold off;

    % Prepare legend text
    legendText = {
        sprintf('l range: [%g, %g] mm', l_min, l_max),
        sprintf('step: %g mm', step),
        sprintf('d: %g mm', d),
        sprintf('\\phi range: [%g, %g] radians', phi_min, phi_max)
    };
    legend(legendText, 'Location', 'bestoutside');
end


%{

function plotResults(positions, phi, kappa, ell)
    figure;
    [X, Y, Z] = cylinder(0.2, 50);
    surf(X, Y, Z, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    hold on;
    scatter3(positions(:,1), positions(:,2), positions(:,3), 'filled');
    
    % Calculate direction vectors
    num_points = size(positions, 1);
    u = zeros(num_points, 1);
    v = zeros(num_points, 1);
    w = zeros(num_points, 1);
    for i = 1:num_points
        theta = kappa(i) * ell(i);
        u(i) = cos(phi(i)) * sin(theta);
        v(i) = sin(phi(i)) * sin(theta);
        w(i) = cos(theta);
    end
    quiver3(positions(:,1), positions(:,2), positions(:,3), u, v, w, 0.5, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 1.5);

    title('Filtered Tip Positions and Directions of the Continuum Robot');
    xlabel('X (mm) - Lateral');
    ylabel('Y (mm) - Longitudinal');
    zlabel('Z (mm) - Vertical');
    grid on;
    axis equal;
    hold off;

     % Prepare legend text
    legendText = {
        sprintf('l range: [%g, %g] mm', l_min, l_max),
        sprintf('d: %g mm', d),
        sprintf('step: %g mm', step),
        sprintf('\\phi range: [%g, %g] radians', phi_min, phi_max)
    };
    legend(legendText, 'Location', 'bestoutside');
end

%}