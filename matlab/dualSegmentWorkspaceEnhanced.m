function dualSegmentWorkspaceEnhanced(r, nb1, nb2, d1, d2, step, t)
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
        last_vector = [cos(last_phi) * sin(last_theta); sin(last_phi) * sin(last_theta); cos(last_theta)];

        % Compute transformation matrix
        T = transformationMatrix(last_phi, last_theta, last_position, last_vector);

        % Create new instance of segment 2
        segment2 = segment2_template;
        for i = 1:size(segment2.positions, 1)
            segment2.positions(i, :) = (T * [segment2.positions(i, :)'; 1]);
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

function T = transformationMatrix(phi, theta, position, vector)
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R = Rz * Ry;
    T = [R, position'; 0 0 0 1];
    translation = [eye(3), vector; 0 0 0 1];
    T = T * translation;
end

function plotNewResults(all_positions, all_phi, all_kappa, all_ell, num_points_segment1)
    figure(1);
    clf;
    hold on;
    scatter3(all_positions(1:num_points_segment1,1), all_positions(1:num_points_segment1,2), all_positions(1:num_points_segment1,3), 36, 'k', 'filled');
    scatter3(all_positions(num_points_segment1+1:end,1), all_positions(num_points_segment1+1:end,2), all_positions(num_points_segment1+1:end,3), 36, all_positions(num_points_segment1+1:end,3), 'filled');
    colormap(jet);
    colorbar;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    title('Enhanced Dual Segments Visualization');
    grid on;
    axis equal;
    hold off;
end
