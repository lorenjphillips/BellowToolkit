function [ss, ms] = ssmslengthcomparison(max_theta, ell)
    
    if mod(max_theta, 2) ~= 0 || max_theta < 0 || max_theta > 180
        error('max_theta must be an even number between 0 and 180.');
    end

    ell_ratio = 50;
    max_theta_ms = max_theta / 2;
    ss = singlesegmentloop(max_theta, ell);
    ms = multisegmentloop(max_theta_ms, ell, ell_ratio);
    compare_ss_ms(ss, ms);

        function compare_ss_ms(ss, ms)
        % Extract angles and distances for both ss and ms
        angles_ss = ss(:, 1);
        distances_ss = ss(:, 2);
        
        angles_ms = ms(:, 1);
        distances_ms = ms(:, 2);
        
        % Create a figure for the scatter plot
        figure;
        hold on;
        
        % Scatter plot for Single Segment (ss)
        scatter(angles_ss, distances_ss, 'b', 'DisplayName', 'Single Segment');
        
        % Scatter plot for Multi Segment (ms)
        scatter(angles_ms, distances_ms, 'r', 'DisplayName', 'Multi Segment');
        
        % Add labels and title
        xlabel('Angle (degrees)');
        ylabel('Distance');
        title('Comparison of Distance at Each Angle for Single and Multi Segment Robots');
        legend;
        grid on;
        
        hold off;
        end

end
