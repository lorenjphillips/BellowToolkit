function g_new = params()
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
    g_new = robotindependentmapping(kappa, phi, ell, 10);
    s = size(g_new);
    s = s(:,1);
    draw_tdcr(g_new,s)
    disp('Result from robot independent mapping:');
    disp(g_new);
end
