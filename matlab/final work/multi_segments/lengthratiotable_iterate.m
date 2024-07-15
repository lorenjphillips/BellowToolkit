function lengthratiotable_iterate(ell)
    % Define ell_ratios to loop through
    ell_ratios = [10, 30, 50, 70, 90];
    curve_ratio = 100;
    max_theta1 = 60;
    max_theta2 = 60;
    max_theta = max_theta1 + max_theta2;
    
    % Compute ss for the given ell value (only once)
    ss = singlesegmentloop(max_theta, ell);
    
    % Initialize a table to store the results
    results_table = table;
    
    for ell_ratio = ell_ratios
        % Compute ms for the given ell value and ell_ratio
        ms = multisegment_iterate(max_theta1, max_theta2, ell, ell_ratio, curve_ratio); % Plots suppressed
        
        % Extract the largest available theta value and the corresponding distances
        [max_theta_value, idx] = max(ms(:, 1)); % Find the index of the largest theta
        distances_at_max_theta = ms(ms(:, 1) == max_theta_value, 2); % Extract distances corresponding to the largest theta
        distance_range = [min(distances_at_max_theta), max(distances_at_max_theta)]; % Determine the range of distances
        max_distance = max(distances_at_max_theta); % Store the maximum distance
        
        % Calculate the average of the ranges of distances for each theta between 61-64º, 65-69º, ..., 116-120º
        theta_ranges = 61:5:116;
        avg_ranges = nan(1, length(theta_ranges));
        
        for k = 1:length(theta_ranges)
            theta_start = theta_ranges(k);
            theta_end = theta_start + 3;
            distances_in_range = ms(ms(:, 1) >= theta_start & ms(:, 1) <= theta_end, 2);
            if ~isempty(distances_in_range)
                avg_range = mean(max(distances_in_range) - min(distances_in_range));
            else
                avg_range = NaN; % Handle case with no distances in range
            end
            avg_ranges(k) = avg_range;
        end
        
        % Append the results to the table
        new_row = [{ell_ratio, max_theta_value, max_distance}, num2cell(avg_ranges)];
        results_table = [results_table; new_row];
    end
    
    % Set the table column names
    col_names = [{'Ell_Ratio', 'Max_Theta', 'Max_Distance'}, strcat('Avg_Range_', string(theta_ranges), '_', string(theta_ranges+3))];
    results_table.Properties.VariableNames = col_names;
    
    % Display the table
    disp(results_table);
    
    % Save the table to a CSV file
    csv_filename = 'results_table.csv';
    writetable(results_table, csv_filename);
    fprintf('Results table saved to %s\n', csv_filename);
    
    % Prompt the user if they want to display the plot
    prompt = 'Do you want to display the plot? Y/N: ';
    user_input = input(prompt, 's');
    
    if strcmpi(user_input, 'Y')
        % Plotting
        figure; fig = figure; fig.Color = [1 1 1];
        hold on;
        
        % Generate colors for different ell ratios
        colors = lines(length(ell_ratios));
        
        % Plot ss results (only once)
        angles_ss = ss(:, 1);
        distances_ss = ss(:, 2);
        scatter(angles_ss, distances_ss, 36, 'k', 'o', 'DisplayName', ['SS - Total Length = ', num2str(ell)]);
        
        for i = 1:length(ell_ratios)
            ell_ratio = ell_ratios(i);
            ms = multisegment_iterate(max_theta1, max_theta2, ell, ell_ratio, curve_ratio); % Plots suppressed
            
            % Plot ms results
            angles_ms = ms(:, 1);
            distances_ms = ms(:, 2);
            scatter(angles_ms, distances_ms, 36, colors(i, :), 'x', 'DisplayName', ['MS Ratio ', num2str(ell_ratio), ':', num2str(100 - ell_ratio)]);
        end
        
        xlabel('Angle (degrees)', 'FontSize', 16);
        ylabel('Distance', 'FontSize', 16);
        title('Distance from Origin at Viewing Angles for Single and Multi Segment', 'FontSize', 14);
        legend('show', 'Location', 'northwest'); % Show legend and set its location
        
        grid on;
        hold off;
    end
end

