function optimize_curves(n, Theta, S_min, S_max, kappa_min, kappa_max)
    % Initial guess for segment lengths and curvatures
    S0 = linspace(S_min/n, S_max/n, n);
    kappa0 = linspace(kappa_min, kappa_max, n);

    % Optimization problem defined
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    x0 = [S0, kappa0]; % Initial guess
    lb = [S_min/n * ones(1, n), kappa_min * ones(1, n)]; % Lower bounds
    ub = [S_max/n * ones(1, n), kappa_max * ones(1, n)]; % Upper bounds

    % Optimization function...
    [x, fval] = fmincon(@(x) -objective(x, n, S_min, S_max), x0, [], [], [], [], lb, ub, @(x) constraints(x, n, Theta, S_min, S_max), options);

    % Extract optimized segment lengths and curvatures
    S_opt = x(1:n);
    kappa_opt = x(n+1:end);

    % Compute curvatures for minimum and maximum lengths
    kappa_min_lengths = kappa_opt .* (S_min ./ sum(S_opt));
    kappa_max_lengths = kappa_opt .* (S_max ./ sum(S_opt));

    % Results
    fprintf('Optimized segment lengths: \n');
    disp(S_opt);
    fprintf('Optimized curvatures: \n');
    disp(kappa_opt);
    fprintf('Curvatures for minimum lengths: \n');
    disp(kappa_min_lengths);
    fprintf('Curvatures for maximum lengths: \n');
    disp(kappa_max_lengths);
    fprintf('Maximum range (X_max - X_min): %f\n', -fval);
end

function f = objective(x, n, S_min, S_max)
    S = x(1:n);
    kappa = x(n+1:end);
    X_max = sum((1 - cos(kappa .* S)) ./ kappa);
    X_min = sum((1 - cos(kappa .* (S / sum(S) * S_min))) ./ kappa);
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
