function maximize_distance_range(S_max, K_max, Theta_target, n, L_min, L_max)
    % Initial guesses for the segment length ratios and curvatures
    S_ratios0 = ones(n, 1) / n;  % Initial guess for arc length ratios
    K0 = linspace(0.01, K_max, n)';  % Initial guess for curvatures, evenly spaced
    x0 = [K0; S_ratios0];

    % Define the objective function
    objective = @(x) -calculate_distance_range(x, n, Theta_target, L_min, L_max);  % Negative because we use fmincon for minimization

    % Define constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [zeros(n, 1); zeros(n, 1)];  % Lower bounds for K and S ratios
    ub = [ones(n, 1) * K_max; ones(n, 1)];  % Upper bounds for K and S ratios (ratio max is 1)

    % Nonlinear constraint to ensure the total length is within bounds and Theta is met
    nonlcon = @(x) nonlinear_constraints(x, Theta_target, n, L_min, L_max);

    % Perform optimization with relaxed tolerances
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                           'StepTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [x_opt, fval, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

    % Extract optimized curvatures and segment length ratios
    K_opt = x_opt(1:n);
    S_ratios_opt = x_opt(n+1:end);

    % Scale the segment lengths to fit within L_min and L_max
    total_length_opt = sum(S_ratios_opt);
    scaling_factor = L_max / total_length_opt;
    S_opt = S_ratios_opt * scaling_factor;

    % Calculate the maximum and minimum distances with the optimized lengths and curvatures
    X_max = calculate_distance(K_opt, S_opt);
    scaling_factor_min = L_min / total_length_opt;
    S_min_opt = S_ratios_opt * scaling_factor_min;
    X_min = calculate_distance(K_opt, S_min_opt);

    % Display results
    disp('Optimized Arc Lengths:');
    disp(S_opt);
    disp('Optimized Curvatures for Max Distance:');
    disp(K_opt);
    disp('Optimized Arc Lengths for Min Distance:');
    disp(S_min_opt);
    disp('Optimized Curvatures for Min Distance:');
    disp(K_opt);
    disp('Maximized Distance Range:');
    disp(X_max - X_min);  % Maximized distance range
    disp('Theta Target (Degrees):');
    disp(rad2deg(Theta_target));  % Convert Theta target to degrees
end

function range = calculate_distance_range(x, n, Theta_target, L_min, L_max)
    K = x(1:n);
    S_ratios = x(n+1:end);

    % Calculate the segment lengths for maximum and minimum distances
    total_length = sum(S_ratios);
    scaling_factor_max = L_max / total_length;
    scaling_factor_min = L_min / total_length;
    S_max = S_ratios * scaling_factor_max;
    S_min = S_ratios * scaling_factor_min;

    % Ensure that the target angle is met
    Theta_max = sum(K .* S_max);
    if abs(Theta_max - Theta_target) > 1e-6
        range = -inf;  % Penalty for not meeting the angle constraint
    else
        X_max = sum((1 - cos(K .* S_max)) ./ K);
        X_min = sum((1 - cos(K .* S_min)) ./ K);
        range = X_max - X_min;
    end
end

function X = calculate_distance(K, S)
    % Calculate distance for given curvatures and arc lengths
    X = sum((1 - cos(K .* S)) ./ K);
end

function [c, ceq] = nonlinear_constraints(x, Theta_target, n, L_min, L_max)
    K = x(1:n);
    S_ratios = x(n+1:end);

    % Calculate the total length
    total_length = sum(S_ratios);
    scaling_factor_max = L_max / total_length;
    S_max = S_ratios * scaling_factor_max;

    % Ensure the target angle is met
    Theta_max = sum(K .* S_max);

    % Inequality constraints for total length
    c = [L_min - total_length; total_length - L_max];

    % Equality constraint to ensure the angle is met
    ceq = Theta_max - Theta_target;
end
% 
% % Example parameters
% S_max = 10;       % Maximum arc length
% K_max = 0.1;      % Maximum curvature
% Theta_target = 1; % Target angle (radians)
% n = 5;            % Number of segments
% L_min = 5;        % Minimum total length
% L_max = 10;       % Maximum total length
% 
% % Run the optimization
% maximize_distance_range(S_max, K_max, Theta_target, n, L_min, L_max);
