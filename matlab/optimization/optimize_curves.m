function optimize_curves(n, Theta, S_min, S_max)
    % Initial guess for segment lengths and curvatures
    % S0_min = linspace(S_min/(n/2), S_min/(2*n), n);
    % S0_max = linspace(S_max/(n/2), S_max/(2*n), n);
    S0_min = [1 ;2 ;3 ;4 ]; % min = 10
    S0_max = [100 ;200 ;300 ;400 ]; % max = 100
    kappa0_min = ones(1, n) * Theta / S_min;
    kappa0_max = ones(1, n) * Theta / S_max;

    % Optimization problem defined
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    x0 = [S0_min, S0_max, kappa0_min, kappa0_max]; % Initial guess
    lb = [0.0 * ones(1, n), 0.0 * ones(1, n), zeros(1, n), zeros(1, n)]; % Lower bounds
    ub = [S_max * ones(1, n), S_max * ones(1, n), inf(1, n), inf(1, n)]; % Upper bounds

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
    ceq3 = sum(S_min_opt) - S_min; % Total segment length for minimum configuration
    ceq4 = sum(S_max_opt) - S_max; % Total segment length for maximum configuration
    ceq = [ceq1; ceq2; ceq3; ceq4];
    
    % Nonlinear inequality constraints
    c1 = kappa_min_opt .* S_min_opt - Theta; % Ensure individual segment curvatures do not exceed viewing angle
    c2 = kappa_max_opt .* S_max_opt - Theta; % Ensure individual segment curvatures do not exceed viewing angle
    
    % Ensure segment lengths stay within reasonable range
    max_segment_length = S_max / 2; % Example upper bound for segment length
    min_segment_length = 0.1; % Example lower bound for segment length
    c3 = S_min_opt - max_segment_length;
    c4 = min_segment_length - S_min_opt;
    c5 = S_max_opt - max_segment_length;
    c6 = min_segment_length - S_max_opt;

    % Ensure all constraints have consistent dimensions
    c = [c1(:); c2(:); c3(:); c4(:); c5(:); c6(:)];
end

% % Example usage
% n = 5;             % Number of segments
% Theta = 1;         % Target viewing angle (radians)
% S_min = 5;         % Minimum total length
% S_max = 10;        % Maximum total length
% 
% optimize_curves(n, Theta, S_min, S_max);
