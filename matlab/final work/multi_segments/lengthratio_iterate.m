function lengthratio_iterate(ell_min, ell_max, ell_steps)
    ell_ratio = 50;
    max_theta1 = 60;
    max_theta2 = 60;
    max_theta = max_theta1 + max_theta1 ;
    curve_ratio = 100;
    ell_values = linspace(ell_min, ell_max, ell_steps);
    
    ss_results = cell(1, ell_steps);
    ms_results = cell(1, ell_steps);
    
    for i = 1:ell_steps % Loop through each ell value and compute ss and ms
        ell = ell_values(i);
        ss = singlesegmentloop(max_theta, ell);
        ms = multisegment_iterate(max_theta1, max_theta2, ell, ell_ratio, curve_ratio); % Plots Supressed
        ss_results{i} = ss; % Save the results in the cell arrays
        ms_results{i} = ms;
    end 
%% Plotting  
figure; fig = figure; fig.Color = [1 1 1];
hold on;
colors = lines(ell_steps); % Generate different colors for each ell value

for i = 1:ell_steps % Loop through the results and plot them
    ss = ss_results{i};
    ms = ms_results{i};
    angles_ss = ss(:, 1);
    distances_ss = ss(:, 2);
    angles_ms = ms(:, 1);
    distances_ms = ms(:, 2);
    
    scatter(angles_ss, distances_ss, 36, colors(i, :), 'o', 'DisplayName', ['SS - Total Length = ', num2str(ell_values(i))]);
    scatter(angles_ms, distances_ms, 36, colors(i, :), 'x', 'DisplayName', ['MS - Total Length = ', num2str(ell_values(i))]);
end

xlabel('Angle (degrees)', 'FontSize', 16);
ylabel('Distance', 'FontSize', 16);
title('Distance from Origin at Viewing Angles for Single and Multi Segment', 'FontSize', 14);
legend('show', 'Location', 'northwest'); % Show legend and set its location
text(0.2, 0.6, ['Segment Length Ratio' ], 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(0.2, 0.55, ['Base:Top = ', num2str(ell_ratio), ':', num2str(abs(ell_ratio - 100))], 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);

grid on;
hold off;


end
