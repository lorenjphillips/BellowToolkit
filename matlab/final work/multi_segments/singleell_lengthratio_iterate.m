function singleell_lengthratio_iterate(ell)
    ell_ratio = 50;
    curve_ratio = 100;
    max_theta1 = 60;
    max_theta2 = 60;
    max_theta = max_theta1 + max_theta2;
    
    % Compute ss and ms for the given ell value
    ss = singlesegmentloop(max_theta, ell);
    ms = multisegment_iterate(max_theta1, max_theta2, ell, ell_ratio, curve_ratio); % Plots suppressed
    
    % Extract the largest available theta value and the corresponding distances
    [max_theta_value, idx] = max(ms(:, 1)); % Find the index of the largest theta
    distances_at_max_theta = ms(ms(:, 1) == max_theta_value, 2); % Extract distances corresponding to the largest theta
    distance_range = [min(distances_at_max_theta), max(distances_at_max_theta)]; % Determine the range of distances
    max_distance = max(distances_at_max_theta); % Store the maximum distance
    
    % Calculate the average of the ranges of distances for each theta between half and full max theta
    half_max_theta = max_theta_value / 2;
    angles_to_consider = ms(:, 1) >= half_max_theta & ms(:, 1) <= max_theta_value;
    unique_thetas = unique(ms(angles_to_consider, 1));
    
    range_sums = 0;
    count = 0;
    
    for theta = unique_thetas'
        distances_at_theta = ms(ms(:, 1) == theta, 2);
        theta_range = max(distances_at_theta) - min(distances_at_theta);
        range_sums = range_sums + theta_range;
        count = count + 1;
    end
    
    avg_distance_range = range_sums / count; % Calculate the average of the distance ranges
    
    % Display the results
    fprintf('Largest available theta: %f degrees\n', max_theta_value);
    fprintf('Range of distances at this theta: [%f, %f]\n', distance_range(1), distance_range(2));
    fprintf('Maximum distance at max theta: %f\n', max_distance);
    fprintf('Average of the ranges of distances between half max theta and max theta: %f\n', avg_distance_range);
    
    % Plotting
    figure; fig = figure; fig.Color = [1 1 1];
    hold on;
    
    % Generate color for the single ell value
    colors = lines(1);
    
    % Plot the results
    angles_ss = ss(:, 1);
    distances_ss = ss(:, 2);
    angles_ms = ms(:, 1);
    distances_ms = ms(:, 2);
    
    scatter(angles_ss, distances_ss, 36, colors(1, :), 'o', 'DisplayName', ['SS - Total Length = ', num2str(ell)]);
    scatter(angles_ms, distances_ms, 36, colors(1, :), 'x', 'DisplayName', ['MS - Total Length = ', num2str(ell)]);
    
    xlabel('Angle (degrees)', 'FontSize', 16);
    ylabel('Distance', 'FontSize', 16);
    title('Distance from Origin at Viewing Angles for Single and Multi Segment', 'FontSize', 14);
    legend('show', 'Location', 'northwest'); % Show legend and set its location
    text(0.2, 0.6, 'Segment Length Ratio', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(0.2, 0.55, ['Base:Top = ', num2str(ell_ratio), ':', num2str(abs(ell_ratio - 100))], 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    grid on;
    hold off;
end
