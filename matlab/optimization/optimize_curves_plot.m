function optimize_curves_plot(n, Theta, S_min, S_max)
    % Initial guess for segment lengths and curvatures
    % Generate random values for min/max initial guesses
    random_min_values = rand(1, n);
    random_max_values = rand(1, n);
    normalized_min_values = random_min_values / sum(random_min_values);
    normalized_max_values = random_max_values / sum(random_max_values);
    S0_min = normalized_min_values * S_min;
    S0_max = normalized_max_values * S_max * 10; % SCALING WILL NOT WORK FOR ALL CASES!!! Will it?

    kappa0_min = ones(1, n) * Theta / S_min;
    kappa0_max = ones(1, n) * Theta / S_max;
        % Adjust bounds if needed
    S_max_bound = S_max/n * 2; % Adjust this as necessary
    S_min_bound = 0.1; % Adjust this as necessary
    
    % Optimization problem defined
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    x0 = [S0_min, S0_max, kappa0_min, kappa0_max]; % Initial guess
    lb = [S_min_bound * ones(1, n), S_min_bound * ones(1, n), zeros(1, n), zeros(1, n)]; % Lower bounds
    ub = [S_max_bound * ones(1, n), S_max_bound * ones(1, n), inf(1, n), inf(1, n)]; % Upper bounds

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

    % Calculate minimum and maximum distances
    X_min = sum((1 - cos(kappa_min .* S_min)) ./ kappa_min);
    X_max = sum((1 - cos(kappa_max .* S_max)) ./ kappa_max);

    % Objective function (maximize X_max - X_min)
    f = -(X_max - X_min);
end

function [c, ceq] = constraints(x, n, Theta, S_min, S_max)
    S_min_opt = x(1:n);
    S_max_opt = x(n+1:2*n);
    kappa_min_opt = x(2*n+1:3*n);
    kappa_max_opt = x(3*n+1:4*n);

    % Nonlinear equality constraints
    ceq1 = sum(kappa_min_opt .* S_min_opt) - Theta; % Maintain viewing angle for minimum configuration
    ceq2 = sum(kappa_max_opt .* S_max_opt) - Theta; % Maintain viewing angle for maximum configuration
    ceq = [ceq1; ceq2];
    
    % Nonlinear inequality constraints
    % Ensure segment lengths stay within reasonable range
    max_segment_length = S_max / 2; % Example upper bound for segment length
    min_segment_length = 0.1; % Example lower bound for segment length
    c3 = S_min_opt - max_segment_length;
    c4 = min_segment_length - S_min_opt;
    c5 = S_max_opt - max_segment_length;
    c6 = min_segment_length - S_max_opt;

    % Ensure all constraints have consistent dimensions
    c = [c3(:); c4(:); c5(:); c6(:)];
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
