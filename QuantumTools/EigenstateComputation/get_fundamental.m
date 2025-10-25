function fundamentalw = get_fundamental(Vf, fieldv, xdomain, Nx, m)

    if nargin < 5
        m = 1;
    end
    Nfield = length(fieldv);
    fundamentalw = zeros(1, Nfield);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
%    V = zeros(Nx, 1);
    % Making a vector of the potential energy at all x:
    if length(Vf) == 1
        % If Vf is a function handle:
        V = Vf(x);
    else
        % If Vf is a vector:
        V = Vf;
    end
%    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
%    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    p = (-pi/dx:(2*pi/xdlength):(pi/dx - 2*pi/xdlength)).';
    % Making a vector of the kinetic energy at all p, for m=1:
    K = p.^2/(2*m);
    Kx = p2x(diag(K));
    for fieldi = 1:Nfield
        Vi = V - fieldv(fieldi)*x;
        H = diag(Vi) + Kx;
        E = eig(H);
        Ereal = sort(real(E));
        fundamentalw(fieldi) = Ereal(2) - Ereal(1);
    end
end