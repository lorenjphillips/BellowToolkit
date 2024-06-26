function visualize_ss_distribution(ss_outputs, angle_bin_size)
    % Combine all angle and distance data from the 3D matrix
    [numSteps, ~, numEll] = size(ss_outputs);
    all_angles = reshape(ss_outputs(:, 1, :), [], 1);
    all_distances = reshape(ss_outputs(:, 2, :), [], 1);

    % Define the angle bins
    angle_min = floor(min(all_angles));
    angle_max = ceil(max(all_angles));
    angle_bins = angle_min:angle_bin_size:angle_max;

    % Initialize the bin counts and distance ranges
    bin_counts = zeros(length(angle_bins)-1, 1);
    distance_min = inf(length(angle_bins)-1, 1);
    distance_max = -inf(length(angle_bins)-1, 1);

    % Bin the angles and calculate the distance ranges
    for i = 1:length(angle_bins)-1
        bin_mask = all_angles >= angle_bins(i) & all_angles < angle_bins(i+1);
        bin_counts(i) = sum(bin_mask);
        if bin_counts(i) > 0
            distance_min(i) = min(all_distances(bin_mask));
            distance_max(i) = max(all_distances(bin_mask));
        end
    end

    % Plot the 3D histogram
    figure;
    Fig = figure;
    Fig.Color = [1 1 1];
    h = bar3(bin_counts);
    for k = 1:length(h)
        zdata = get(h(k),'ZData');
        set(h(k),'CData',zdata,'FaceColor','interp')
    end
    colorbar;
    xlabel('Angle Bins (degrees)');
    ylabel('Count');
    zlabel('Number of Occurrences');
    set(gca, 'XTick', 1:length(angle_bins)-1, 'XTickLabel', angle_bins(1:end-1));
    title('Distribution of Angles and Distances');
    
    % Adjust color scale to max occurrences
    max_occurrences = max(bin_counts);
    clim([0 max_occurrences]);

    % Prepare data for the table
    table_data = cell(length(angle_bins)-1, 2);
    for i = 1:length(angle_bins)-1
        if bin_counts(i) > 0
            table_data{i, 1} = sprintf('%.2f - %.2f', distance_min(i), distance_max(i));
            table_data{i, 2} = sprintf('%.2f', angle_bins(i));
        else
            table_data{i, 1} = '';
            table_data{i, 2} = '';
        end
    end

    % Create table annotation
    annotation('textbox', [0.7, 0.5, 0.25, 0.4], 'String', ...
        {'Distance Range', 'Angle (degrees)', table_data{:,:}}, ...
        'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'FitBoxToText', 'off');

end
