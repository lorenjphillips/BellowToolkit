function lengthratio_iterate(ell_min, ell_max, ell_steps)
    % if mod(max_theta, 2) ~= 0 || max_theta < 0 || max_theta > 180
    %     error('max_theta must be an even number between 0 and 180.');
    % end
    
    ell_ratio = 75;
    max_theta1 = 50;
    max_theta2 = 50;
    max_theta = max_theta1 + max_theta1 ;
    curve_ratio = 100;
    
    % Generate the range of ell values
    ell_values = linspace(ell_min, ell_max, ell_steps);
    
    % Initialize cell arrays to store results
    ss_results = cell(1, ell_steps);
    ms_results = cell(1, ell_steps);
    
    % Loop through each ell value and compute ss and ms
    for i = 1:ell_steps
        ell = ell_values(i);
        ss = singlesegmentloop(max_theta, ell);
        % ms = multisegmentloop_ratio(max_theta_ms, ell, ell_ratio);
        ms = multisegmentloop_ratio_noplots(max_theta1, max_theta2, ell, ell_ratio, curve_ratio);
        %SUPRESS THE PLOTS 
        % Save the results in the cell arrays
        ss_results{i} = ss;
        ms_results{i} = ms;
    end
  %% Plotting  
    % Plot the results
    figure; fig = figure; fig.Color = [1 1 1];
    hold on;
    
    colors = lines(ell_steps); % Generate different colors for each ell value
    
    % Loop through the results and plot them
    for i = 1:ell_steps
        ss = ss_results{i};
        ms = ms_results{i};
        
        angles_ss = ss(:, 1);
        distances_ss = ss(:, 2);
        
        angles_ms = ms(:, 1);
        distances_ms = ms(:, 2);
        
        scatter(angles_ss, distances_ss, 36, colors(i, :), 'o', 'DisplayName', ['SS - ell = ', num2str(ell_values(i))]);
        scatter(angles_ms, distances_ms, 36, colors(i, :), 'x', 'DisplayName', ['MS - ell = ', num2str(ell_values(i))]);
    end
    
    xlabel('Angle (degrees)');
    ylabel('Distance');
    title('Comparison of Distance at Each Angle for Single and Multi Segments');
    legend show;
    grid on;
    
    hold off;
end
