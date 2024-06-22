function [phi_val, ell_val, kappa_val] = continuumRobotKinematics(l1, l2, l3, d)

    % Calculate the orientation angle phi(q)
    phi_val = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));

    % Calculate the effective length ell(q)
    ell_val = (l1 + l2 + l3) / 3;

    % Calculate the curvature kappa(q)
    kappa_val = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));

end