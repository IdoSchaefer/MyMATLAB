function performance = percosVsym(Vk, xdomain, kdomain, Nk, penalty)
    sz = size(Vk);
    if sz(1) == 1
        Ncosines = sz(2)/2;
        Vk = Vk.';
    else
        Ncosines = sz(1)/2;
    end
    Nx = 2*Ncosines - 1;
    L = xdomain(2) - xdomain(1);
    dx = L/(Ncosines - 1);
    Vx = zeros(Nx + 2, 1);
    % The first and last term will remain 0.
    Vx(2:((Nx + 1)/2 + 1)) = dctI(Vk(1:Ncosines)) + 1i*dctI(Vk((Ncosines + 1):2*Ncosines));
    Vx(((Nx + 1)/2 + 2):(Nx + 1)) = Vx(((Nx + 1)/2):-1:2);
    %performance = Vabs_efficiency(Vx, [xdomain(1) - dx, xdomain(2) + L + dx], kdomain, Nk) + penalty*sum(Vk.^2);
    performance = Vabs_efficiency1(Vx, [xdomain(1), xdomain(2) + L], kdomain, Nk) + penalty*sum(Vk.^2);
end        