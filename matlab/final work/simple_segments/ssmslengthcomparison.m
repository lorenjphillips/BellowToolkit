function [ss, ms] = ssmslengthcomparison(max_theta, ell)
    
    if mod(max_theta, 2) ~= 0 || max_theta < 0 || max_theta > 180
        error('max_theta must be an even number between 0 and 180.');
    end

    ell_ratio = 50;
    max_theta_ms = max_theta / 2;
    ss = singlesegmentloop(max_theta, ell);
    ms = multisegmentloop(max_theta_ms, ell, ell_ratio);
end
