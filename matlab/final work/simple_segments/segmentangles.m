function  segmentangles()
    n = input('Enter the number of segments: ');
    kappa = zeros(1, n);
    phi = zeros(1, n);
    ell = zeros(1, n);
    for i = 1:n
        fprintf('Enter data for segment %d:\n', i);
        l1 = input('Enter l1: ');
        l2 = input('Enter l2: ');
        l3 = input('Enter l3: ');
        d = input('Enter d: ');
        kappa(i) = 2 * sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3) / (d * (l1 + l2 + l3));
        phi(i) = atan2(sqrt(3) * (l2 + l3 - 2 * l1), 3 * (l2 - l3));
        ell(i) = (l1 + l2 + l3) / 3;
    end
g = robotindependentmapping(kappa, phi, ell, 20);
assignin('base','g',g)
    
% Extract the components of the vector from the last row of g
vx = g(end, 9);
vy = g(end, 10);
vz = g(end, 11);

% Calculate the magnitude of the vector
v_mag = sqrt(vx^2 + vy^2 + vz^2);

% Calculate the magnitude of the projection of the vector onto the x-y plane
v_proj_mag = sqrt(vx^2 + vy^2);

% Calculate the angle between the vector and the x-y plane
angle_rad = acos(v_proj_mag / v_mag);

% Convert the angle from radians to degrees
angle_deg = rad2deg(angle_rad);

% Display the result
disp(['The angle with respect to the x-y plane is ', num2str(angle_deg), ' degrees']);
disp(g)
end
