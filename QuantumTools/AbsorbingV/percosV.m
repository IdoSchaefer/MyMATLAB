function performance = percosV(Vk, xdomain, kdomain, Nk, penalty)
% The function computes the performance parameter for a cosine series
% absorbing potential.
    sz = size(Vk);
    if sz(1) == 1
        Nx = sz(2)/2;
        Vk = Vk.';
    else
        Nx = sz(1)/2;
    end
    %dx = (xdomain(2) - xdomain(1))/(Nx - 1);
    %Vx = [0; dctI(Vk(1:Nx)) + 1i*dctI(Vk((Nx + 1):2*Nx)); 0];
    Vx = dctI(Vk(1:Nx)) + 1i*dctI(Vk((Nx + 1):2*Nx));
%    performance = Vabs_efficiency(Vx, [xdomain(1) - dx, xdomain(2) + dx], kdomain, Nk) + penalty*sum(Vk.^2);
    performance = Vabs_efficiency1(Vx, xdomain, kdomain, Nk) + penalty*sum(Vk.^2);
end