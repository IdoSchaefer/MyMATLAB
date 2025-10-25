function viewPp(U, dx, xdlength, dt)
% The function views the evolution of the wave function in the p domain.
    if nargin<3
        dt = 0.1;
    end
    szU = size(U);
    Nx = szU(1);
    Nt = szU(2);
    p = (-pi/dx:(2*pi/xdlength):(pi*(1/dx - 1/xdlength))).';
    Up = fft(U)/sqrt(Nx);
    Up = [Up((Nx/2 + 1):Nx, :); Up(1:Nx/2, :)];
    Pp = conj(Up).*Up;
    figure
    axis([p(1), p(Nx), min(min(Pp)), max(max(Pp))])
    Pcurve = line(p, Pp(:, 1));
    for ti = 2:Nt
        pause(dt);
        set(Pcurve, 'Ydata', Pp(:, ti))
    end
end