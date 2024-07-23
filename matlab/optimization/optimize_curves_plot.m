function optimize_curves_plot(n, Theta, S_min, S_max, kappa_min, kappa_max)
    % Initial guess for segment lengths and curvatures
    S0 = linspace(S_min/n, S_max/n, n);
    kappa0 = linspace(kappa_min, kappa_max, n);

    % Optimization problem defined
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    x0 = [S0, kappa0]; % Initial guess
    lb = [S_min/n * ones(1, n), kappa_min * ones(1, n)]; % Lower bounds
    ub = [S_max/n * ones(1, n), kappa_max * ones(1, n)]; % Upper bounds

    % Optimization function
    [x, fval] = fmincon(@(x) -objective(x, n, S_min, S_max), x0, [], [], [], [], lb, ub, @(x) constraints(x, n, Theta, S_min, S_max), options);

    % Extract optimized segment lengths and curvatures
    S_opt = x(1:n);
    kappa_opt = x(n+1:end);

    % Compute minimum and maximum lengths for each segment
    scaling_factor_max = S_max / sum(S_opt);
    scaling_factor_min = S_min / sum(S_opt);
    S_max_lengths = S_opt * scaling_factor_max;
    S_min_lengths = S_opt * scaling_factor_min;

    % Results
    fprintf('Optimized segment lengths: \n');
    disp(S_opt);
    fprintf('Optimized curvatures: \n');
    disp(kappa_opt);
    fprintf('Segment lengths for minimum total length: \n');
    disp(S_min_lengths);
    fprintf('Segment lengths for maximum total length: \n');
    disp(S_max_lengths);
    fprintf('Maximum range (X_max - X_min): %f\n', -fval);

    % Compute g_min and g_max using robotindependentmapping function
    phi = zeros(n, 1); % phi is 0 for all cases
    ptsperseg = 20; % Number of points per segment

    g_min = robotindependentmapping(kappa_opt, phi, S_min_lengths, ptsperseg);
    g_max = robotindependentmapping(kappa_opt, phi, S_max_lengths, ptsperseg);

    % Display results
    fprintf('g_min: \n');
    disp(g_min);
    fprintf('g_max: \n');
    disp(g_max);
end

function f = objective(x, n, S_min, S_max)
    S = x(1:n);
    kappa = x(n+1:end);
    
    % Calculate maximum distance
    scaling_factor_max = S_max / sum(S);
    S_max_lengths = S * scaling_factor_max;
    X_max = sum((1 - cos(kappa .* S_max_lengths)) ./ kappa);
    
    % Calculate minimum distance
    scaling_factor_min = S_min / sum(S);
    S_min_lengths = S * scaling_factor_min;
    X_min = sum((1 - cos(kappa .* S_min_lengths)) ./ kappa);
    
    f = X_max - X_min; % Objective function (maximize X_max - X_min)
end

function [c, ceq] = constraints(x, n, Theta, S_min, S_max)
    S = x(1:n);
    kappa = x(n+1:end);
    
    % Nonlinear equality constraints
    ceq = sum(kappa .* S) - Theta; % Maintain viewing angle
    
    % Nonlinear inequality constraints
    c = []; % None in this case
end

% % Example usage
% n = 5;             % Number of segments
% Theta = 1;         % Target viewing angle (radians)
% S_min = 5;         % Minimum total length
% S_max = 10;        % Maximum total length
% kappa_min = 0.01;  % Minimum curvature
% kappa_max = 0.1;   % Maximum curvature
% 
% optimize_curves_plot(n, Theta, S_min, S_max, kappa_min, kappa_max);

