function continuumKinematics(l1_min, l1_max, d, step)
    % Validate input parameters
    if l1_min >= l1_max
        error('Minimum tendon length must be less than maximum tendon length');
    end
    if step <= 0
        error('Step must be a positive value');
    end
    if d <= 0
        error('Diameter must be a positive value');
    end

    % Initialize the range for each tendon
    l_range = l1_min:step:l1_max;
    num_elements = numel(l_range);

    % Preallocate arrays for calculated values
    phi = zeros(num_elements^3, 1);
    ell = zeros(num_elements^3, 1);
    kappa = zeros(num_elements^3, 1);
    lengths = zeros(num_elements^3, 3);
    positions = zeros(num_elements^3, 3);

    % Iterate over each tendon length
    index = 1;
    for l1 = l_range
        for l2 = l_range
            for l3 = l_range
                % Kinematics calculations
                [current_phi, current_ell, current_kappa] = calculateKinematics(l1, l2, l3, d);

                % Append results
                phi(index) = current_phi;
                ell(index) = current_ell;
                kappa(index) = current_kappa;
                lengths(index, :) = [l1, l2, l3];

                % Position calculations
                positions(index, :) = calculatePosition(current_ell, current_kappa, current_phi);

                index = index + 1;
            end
        end
    end

    % Plotting results
    plotResults(positions);
end

function [phi, ell, kappa] = calculateKinematics(l1, l2, l3, d)
    phi = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));
    ell = (l1 + l2 + l3) / 3;
    kappa = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));
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
        z = ell; % No bending, only extension
        x = 0;
        y = 0;
    end
    position = [x, y, z];
end

function plotResults(positions)
    figure(1);
    f = figure;
    f.Color = [1 1 1];              % Set the color to white
    [x, y, z] = cylinder(0.2, 5);
    surf(x, y, z, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    hold on;
    scatter3(positions(:,1), positions(:,2), positions(:,3), 36, positions(:,3), 'filled');
    text(0, 0, 1.1, 'Base', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    title('Tip Positions of the Continuum Robot', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X (mm) - Lateral', 'FontSize', 12);
    ylabel('Y (mm) - Longitudinal', 'FontSize', 12);
    zlabel('Z (mm) - Along Bending Direction', 'FontSize', 12);
    colormap(jet);
    colorbar('Label', 'Height (Z) in mm');
    grid on;
    axis equal;
    hold off;
end