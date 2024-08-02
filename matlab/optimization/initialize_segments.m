function [S0_min, S0_max] = initialize_segments(n, S_min, S_max, ratio_type, percentage)
    if nargin < 5
        percentage = 0; % Default percentage if not provided
    end

    switch ratio_type
        case 'increasing'
            % Calculate the segment lengths
            increments = (1 + percentage / 100) .^ (0:n-1);
            S0_min_ratios = increments / sum(increments);
        case 'decreasing'
            % Calculate the segment lengths
            increments = (1 - percentage / 100) .^ (0:n-1);
            S0_min_ratios = increments / sum(increments);
        case 'equal'
            % Equal lengths
            S0_min_ratios = ones(1, n) / n;
        otherwise
            error('Invalid ratio_type. Use "increasing", "decreasing", or "equal".');
    end

    % Set the segment lengths according to the ratios
    S0_min = S0_min_ratios * S_min;
    S0_max = S0_min_ratios * S_max;
end
