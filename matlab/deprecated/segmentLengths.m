function [l_1_min, l_1_max, l_2_min, l_2_max] = segmentLengths(r, nb1, nb2, t)
    % t = 1.5; % Thickness in mm
    l_1_min = t * nb1;
    l_1_max = nb1 * 2 * pi * r;
    l_2_min = t * nb2;
    l_2_max = nb2 * 2 * pi * r;
end