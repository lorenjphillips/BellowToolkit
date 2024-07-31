function optimize_curves_plot_cs(n, Theta, S_min, S_max)
    % Initial guess for segment lengths and curvatures
    % Generate random values for min/max initial guesses
    random_min_values = rand(1, n);
    random_max_values = rand(1, n);
    normalized_min_values = random_min_values / sum(random_min_values);
    normalized_max_values = random_max_values / sum(random_max_values);
    S0_min = normalized_min_values * S_min;
    S0_max = normalized_max_values * S_max;

    % Initial guesses for curvatures
    kappa0_min = max(rand(1, n) * Theta / S_min, 1e-6); % Avoid zero curvatures
    ratio_kappa = kappa0_min / kappa0_min(1); % Ratio of curvatures

    % Compute corresponding maximum curvatures maintaining the ratio
    kappa0_max = ratio_kappa * (Theta / S_max) * (S_max / S_min);

    % Optimization problem defined
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                           'MaxFunctionEvaluations', 10000, 'ConstraintTolerance', 1e-6);

    % Ensure x0 is within the bounds
    x0 = [S0_min, S0_max, kappa0_min, kappa0_max]; % Initial guess

    % Define lower and upper bounds
    lb = [zeros(1, n), zeros(1, n), zeros(1, n), zeros(1, n)]; % Lower bounds
    ub = [ones(1, n) * S_min, ones(1, n) * S_max, inf(1, n), inf(1, n)]; % Upper bounds

    % Adjust x0 to fit within lb and ub
    x0 = max(lb, min(x0, ub));

    % Optimization function
    [x, fval] = fmincon(@objective, x0, [], [], [], [], lb, ub, @(x) constraints(x, n, Theta, S_min, S_max), options);

    % Extract optimized segment lengths and curvatures for minimum and maximum configurations
    S_min_opt = x(1:n);
    S_max_opt = x(n+1:2*n);
    kappa_min_opt = x(2*n+1:3*n);
    kappa_max_opt = x(3*n+1:4*n);

    % Compute minimum and maximum distances
    X_min = sum((1 - cos(kappa_min_opt .* S_min_opt)) ./ kappa_min_opt);
    X_max = sum((1 - cos(kappa_max_opt .* S_max_opt)) ./ kappa_max_opt);

    % Results
    fprintf('Optimized segment lengths for minimum configuration: \n');
    disp(S_min_opt);
    fprintf('Optimized curvatures for minimum configuration: \n');
    disp(kappa_min_opt);
    fprintf('Optimized segment lengths for maximum configuration: \n');
    disp(S_max_opt);
    fprintf('Optimized curvatures for maximum configuration: \n');
    disp(kappa_max_opt);
    fprintf('Maximum range (X_max - X_min): %f\n', X_max - X_min);

    % Compute g_min and g_max using robotindependentmapping function
    phi = zeros(n, 1); % phi is 0 for all cases
    ptsperseg = 20; % Number of points per segment

    g_min = robotindependentmapping(kappa_min_opt, phi, S_min_opt, ptsperseg);
    g_max = robotindependentmapping(kappa_max_opt, phi, S_max_opt, ptsperseg);

    % Plot the results
    plot_robot_segments(g_min, g_max, n);
end

function f = objective(x)
    n = length(x) / 4;
    S_min = x(1:n);
    S_max = x(n+1:2*n);
    kappa_min = x(2*n+1:3*n);
    kappa_max = x(3*n+1:4*n);

    % Avoid division by zero or near-zero values
    kappa_min = max(kappa_min, 1e-6);
    kappa_max = max(kappa_max, 1e-6);

    % Calculate minimum and maximum distances
    X_min = sum((1 - cos(kappa_min .* S_min)) ./ kappa_min);
    X_max = sum((1 - cos(kappa_max .* S_max)) ./ kappa_max);

    % Objective function (maximize X_max - X_min)
    f = -(X_max - X_min)^2;
end

function [c, ceq] = constraints(x, n, Theta, S_min, S_max)
    S_min_opt = x(1:n);
    S_max_opt = x(n+1:2*n);
    kappa_min_opt = x(2*n+1:3*n);
    kappa_max_opt = x(3*n+1:4*n);

    % Avoid division by zero or near-zero values
    kappa_min_opt = max(kappa_min_opt, 1e-6);
    kappa_max_opt = max(kappa_max_opt, 1e-6);

    % Nonlinear equality constraints
    ceq1 = sum(kappa_min_opt .* S_min_opt) - Theta; % Maintain viewing angle for minimum configuration
    ceq2 = sum(kappa_max_opt .* S_max_opt) - Theta; % Maintain viewing angle for maximum configuration
    ceq3 = kappa_max_opt / kappa_max_opt(1) - kappa_min_opt / kappa_min_opt(1); % Maintain curvature ratio
    ceq4 = S_max_opt / S_max_opt(1) - S_min_opt / S_min_opt(1); % Maintain length ratio
    ceq5 = sum(S_max_opt) - S_max; % Total length of maximum configuration must be S_max
    ceq = [ceq1; ceq2; ceq3(:); ceq4(:); ceq5];
    
    % Nonlinear inequality constraints
    c1 = 1 ./ kappa_min_opt - S_min_opt / pi; % Limit curvature for minimum configuration
    c2 = 1 ./ kappa_max_opt - S_max_opt / pi; % Limit curvature for maximum configuration
    c3 = S_min - sum(S_min_opt); % Total length of minimum configuration must be >= S_min

    % Ensure all constraints have consistent dimensions
    c = [c1(:); c2(:); c3];
end

function plot_robot_segments(g_min, g_max, n)
    figure;
    fig = figure;
    fig.Color = [1 1 1];
    hold on;
    colors = lines(n); % Generate n distinct colors

    % Plot minimum lengths configuration
    for idx = 1:n
        seg_end = idx * 20; % Assuming ptsperseg = 20
        if idx == 1
            seg_start = 1;
        else
            seg_start = (idx - 1) * 20 + 1;
        end
        
        % Extract the origin and components of the vector
        x_min = g_min(seg_start:seg_end, 13); 
        z_min = g_min(seg_start:seg_end, 15);
        
        % Plot the backbone
        plot(x_min, z_min, 'LineWidth', 2, 'Color', colors(idx, :));
    end

    % Plot maximum lengths configuration
    for idx = 1:n
        seg_end = idx * 20; % Assuming ptsperseg = 20
        if idx == 1
            seg_start = 1;
        else
            seg_start = (idx - 1) * 20 + 1;
        end
        
        % Extract the origin and components of the vector
        x_max = g_max(seg_start:seg_end, 13); 
        z_max = g_max(seg_start:seg_end, 15);
        
        % Plot the backbone
        plot(x_max, z_max, 'LineWidth', 2, 'Color', colors(idx, :));
    end

    xlabel('X');
    ylabel('Y');
    title('Plot of Minimum and Maximum Lengths');
    grid on;
    hold off;
end

% % Example usage
% n = 5;             % Number of segments
% Theta = pi;        % Target viewing angle (radians)
% S_min = 5;         % Minimum total length
% S_max = 10;        % Maximum total length
% 
% optimize_curves_plot_cs(n, Theta, S_min, S_max);
