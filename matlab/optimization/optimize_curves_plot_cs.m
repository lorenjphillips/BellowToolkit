function optimize_curves_plot_cs(n, Theta, S_min, S_max, ratio_type, percentage)
    % Initialize fixed segment lengths
    [S_min_fixed, S_max_fixed] = initialize_segments(n, S_min, S_max, ratio_type, percentage);

    % Initial guess for curvatures
    kappa0_min = initialize_curvatures(S_min_fixed, Theta);
    kappa0_max = initialize_curvatures(S_max_fixed, Theta);


    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', 10000, 'StepTolerance', 1e-9, ...
                       'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
                       'OutputFcn', @printCurvatures);

    function stop = printCurvatures(x,optimValues,state)
    stop = false;
    n = length(x) / 2;
    kappa_min = x(1:n);
    kappa_max = x(n+1:end);

    % Print curvatures at each iteration
    fprintf('Iteration %d: Min Curvatures = %s, Max Curvatures = %s\n', ...
            optimValues.iteration, mat2str(kappa_min), mat2str(kappa_max));
    end


    % Initial guess combining curvatures only
    x0 = [kappa0_min, kappa0_max]; 

    % Define lower and upper bounds for curvatures
    lb = zeros(1, 2*n); % Lower bounds
    ub = inf(1, 2*n);   % Upper bounds

    % Adjust x0 to fit within lb and ub
    x0 = max(lb, min(x0, ub));

    % Optimization function
    [x, fval] = fmincon(@(x) objective(x, S_min_fixed, S_max_fixed), x0, [], [], [], [], lb, ub, @(x) constraints(x, n, Theta, S_min_fixed, S_max_fixed), options);

    % Extract optimized curvatures for minimum and maximum configurations
    kappa_min_opt = x(1:n);
    kappa_max_opt = x(n+1:2*n);

    % Compute minimum and maximum distances
    X_min = sum((1 - cos(kappa_min_opt .* S_min_fixed)) ./ kappa_min_opt);
    X_max = sum((1 - cos(kappa_max_opt .* S_max_fixed)) ./ kappa_max_opt);

    % Results
    fprintf('Fixed segment lengths for minimum configuration: \n');
    disp(S_min_fixed);
    fprintf('Optimized curvatures for minimum configuration: \n');
    disp(kappa_min_opt);
    fprintf('Fixed segment lengths for maximum configuration: \n');
    disp(S_max_fixed);
    fprintf('Optimized curvatures for maximum configuration: \n');
    disp(kappa_max_opt);
    fprintf('Maximum range (X_max - X_min): %f\n', X_max - X_min);

    % Compute g_min and g_max using robotindependentmapping function
    phi = zeros(n, 1); % phi is 0 for all cases
    ptsperseg = 20; % Number of points per segment

    g_min = robotindependentmapping(kappa_min_opt, phi, S_min_fixed, ptsperseg);
    g_max = robotindependentmapping(kappa_max_opt, phi, S_max_fixed, ptsperseg);

    % Plot the results along with initial guesses
    plot_robot_segments(g_min, g_max, n, kappa0_min, kappa0_max, S_min_fixed, S_max_fixed);
end

function f = objective(x, S_min_fixed, S_max_fixed)
    n = length(x) / 2;
    kappa_min = x(1:n);
    kappa_max = x(n+1:2*n);

    % Avoid division by zero or near-zero values
    kappa_min = max(kappa_min, 1e-6);
    kappa_max = max(kappa_max, 1e-6);

    % Calculate minimum and maximum distances
    X_min = sum((1 - cos(kappa_min .* S_min_fixed)) ./ kappa_min);
    X_max = sum((1 - cos(kappa_max .* S_max_fixed)) ./ kappa_max);

    % Objective function (maximize X_max - X_min)
    f = -(X_max - X_min)^2;
end

function kappa0 = initialize_curvatures(S, Theta)
    n = length(S);

    % Each curvature should be equal to the target angle divided by the number of segments
    kappa0 = repmat(Theta / n, 1, n);

    % Ensure curvatures are within physical limits (i.e., bending radius constraints)
    max_curvature = pi ./ S;
    kappa0 = min(kappa0, max_curvature);
end


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

function [c, ceq] = constraints(x, n, Theta, S_min_fixed, S_max_fixed)
    kappa_min_opt = x(1:n);
    kappa_max_opt = x(n+1:2*n);

    % Ensure curvatures are positive
    kappa_min_opt = max(kappa_min_opt, 1e-6);
    kappa_max_opt = max(kappa_max_opt, 1e-6);

    % Nonlinear equality constraints
    ceq1 = sum(kappa_min_opt .* S_min_fixed) - Theta; % Ensure total curvature for min configuration equals Theta
    ceq2 = sum(kappa_max_opt .* S_max_fixed) - Theta; % Ensure total curvature for max configuration equals Theta

    % Nonlinear inequality constraints (unchanged)
    c1 = 1 ./ kappa_min_opt - S_min_fixed / pi;
    c2 = 1 ./ kappa_max_opt - S_max_fixed / pi;
    c3 = sum(kappa_min_opt .* S_min_fixed) - pi; % Ensure no over-bending
    c4 = sum(kappa_max_opt .* S_max_fixed) - pi; % Ensure no over-bending

    ceq = [ceq1; ceq2];
    c = [c1(:); c2(:); c3; c4];
end


function plot_robot_segments(g_min, g_max, n, kappa0_min, kappa0_max, S_min_fixed, S_max_fixed)
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

    % Plot initial guesses for minimum and maximum configurations
    phi = zeros(n, 1); % phi is 0 for all cases
    ptsperseg = 20; % Number of points per segment
    g_min_initial = robotindependentmapping(kappa0_min, phi, S_min_fixed, ptsperseg);
    g_max_initial = robotindependentmapping(kappa0_max, phi, S_max_fixed, ptsperseg);

    % Plot initial guess for minimum configuration
    for idx = 1:n
        seg_end = idx * 20;
        if idx == 1
            seg_start = 1;
        else
            seg_start = (idx - 1) * 20 + 1;
        end

        % Extract the origin and components of the vector
        x_min_init = g_min_initial(seg_start:seg_end, 13); 
        z_min_init = g_min_initial(seg_start:seg_end, 15);

        % Plot the initial guess backbone
        plot(x_min_init, z_min_init, '--', 'LineWidth', 1, 'Color', colors(idx, :));
        text(x_min_init(end), z_min_init(end), sprintf('Min Init %d', idx), 'FontSize', 8);
    end

    % Plot initial guess for maximum configuration
    for idx = 1:n
        seg_end = idx * 20;
        if idx == 1
            seg_start = 1;
        else
            seg_start = (idx - 1) * 20 + 1;
        end

        % Extract the origin and components of the vector
        x_max_init = g_max_initial(seg_start:seg_end, 13); 
        z_max_init = g_max_initial(seg_start:seg_end, 15);

        % Plot the initial guess backbone
        plot(x_max_init, z_max_init, '--', 'LineWidth', 1, 'Color', colors(idx, :));
        text(x_max_init(end), z_max_init(end), sprintf('Max Init %d', idx), 'FontSize', 8);
    end

    xlabel('X');
    ylabel('Y');
    title('Plot of Minimum and Maximum Lengths with Initial Guesses');
    grid on;
    hold off;
end
