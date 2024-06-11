function position = calculatePosition(ell, kappa, phi)
    if kappa ~= 0
        R = 1 / kappa; % Radius of curvature
        theta = kappa * ell; % Bend angle
        dist = R - (R * cos(theta)); % Distance from origin
        z = R * sin(theta); % along the bending direction
        x = dist * cos(phi); % lateral displacement
        y = dist * sin(phi); % depth displacement
    else
        x = ell * cos(phi); % If no curvature, position extends in straight line
        y = ell * sin(phi);
        z = 0; % No bending, hence Z remains constant
    end
    position = [x, y, z];
end
