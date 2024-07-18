function maximize_distance_range(S_max, K_max, Theta_target, n, L_min, L_max)
    % Initial guesses for the arc lengths
    S0 = linspace(0.01, S_max, n)';  % Initial guess for arc lengths, evenly spaced

    % Define the objective function for arc length optimization
    objective_S = @(S) -distance_range_with_fixed_theta(S, K_max, Theta_target, n);  % Negative because we use fmincon for minimization

    % Define constraints for arc length optimization
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = zeros(n, 1);  % Lower bounds for S
    ub = ones(n, 1) * S_max;  % Upper bounds for S
    nonlcon_S = @(S) nonlinear_constraints_length(S, L_min, L_max);

    % Perform optimization of arc lengths with relaxed tolerances
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                           'StepTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [S_opt, fval, exitflag, output] = fmincon(objective_S, S0, A, b, Aeq, beq, lb, ub, nonlcon_S, options);

    % Find the optimized curvatures for max and min distances
    [X_max, X_min, K_max_opt, K_min_opt] = distance_range_with_fixed_theta(S_opt, K_max, Theta_target, n);

    % Display results
    disp('Optimized Arc Lengths:');
    disp(S_opt);
    disp('Optimized Curvatures for Max Distance:');
    disp(K_max_opt);
    disp('Optimized Curvatures for Min Distance:');
    disp(K_min_opt);
    disp('Maximized Distance Range:');
    disp(X_max - X_min);  % Maximized distance range
    disp('Theta Target (Degrees):');
    disp(rad2deg(Theta_target));  % Convert Theta target to degrees
end

function [X_max, X_min, K_max_opt, K_min_opt] = distance_range_with_fixed_theta(S, K_max, Theta_target, n)
    % Calculate curvatures that achieve the target Theta with given S
    K_max_opt = Theta_target ./ S;
    K_min_opt = zeros(n, 1);  % Minimum curvature is zero

    % Ensure the curvatures do not exceed K_max
    K_max_opt = min(K_max_opt, K_max);

    % Calculate maximum and minimum distances
    X_max = sum((1 - cos(K_max_opt .* S)) ./ K_max_opt);
    X_min = sum((1 - cos(K_min_opt .* S)) ./ K_min_opt);  % Should be zero

end

function [c, ceq] = nonlinear_constraints_length(S, L_min, L_max)
    % Calculate the total length
    total_length = sum(S);

    % Inequality constraints for total length
    c = [L_min - total_length; total_length - L_max];

    % No equality constraints
    ceq = [];
end


%% Example parameters
% S_max = 10;       % Maximum arc length
% K_max = 0.1;      % Maximum curvature
% Theta_target = 1; % Target angle (radians)
% n = 5;            % Number of segments
% L_min = 5;        % Minimum total length
% L_max = 10;       % Maximum total length
% 
% % Run the optimization
% maximize_distance_range(S_max, K_max, Theta_target, n, L_min, L_max);
