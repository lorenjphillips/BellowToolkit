function segment = calculateSegment(l_min, l_max, d, step)
    l_range = l_min:step:l_max;
    num_elements = numel(l_range);
    phi = zeros(num_elements^3, 1);
    ell = zeros(num_elements^3, 1);
    kappa = zeros(num_elements^3, 1);
    positions = zeros(num_elements^3, 3);
    index = 1;
    for l1 = l_range
        for l2 = l_range
            for l3 = l_range
                [phi_val, ell_val, kappa_val] = continuumRobotKinematics(l1, l2, l3, d);
                %[phi_val, ell_val, kappa_val] = calculateKinematics(l1, l2, l3, d);
                phi(index) = phi_val;
                ell(index) = ell_val;
                kappa(index) = kappa_val;
                positions(index, :) = calculatePosition(ell_val, kappa_val, phi_val);
                index = index + 1;
            end
        end
    end
    segment = struct('positions', positions, 'phi', phi, 'kappa', kappa, 'ell', ell);
end