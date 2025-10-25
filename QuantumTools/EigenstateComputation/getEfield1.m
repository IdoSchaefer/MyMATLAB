function [Efield, es_field]  = getEfield1(Vf, miu, fieldv, n_levels, xdomain, m)

    if nargin < 6
        m = 1;
    end
    Nfield = length(fieldv);
    Nx = length(miu);
    Nlevels = length(n_levels);
    Efield = zeros(Nlevels, Nfield);
    es_field = zeros(Nx, Nfield, Nlevels);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
%    V = zeros(Nx, 1);
    % Creating a vector of the potential energy at all x:
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
    Hfieldi = diag(V) + Kx;
    [Pfieldi, D] = eig(Hfieldi);
    Efieldi = diag(D);
    [Efieldi, orderE] = sort(real(Efieldi));
    Pfieldi = Pfieldi(:, orderE);
    Efield(:, 1) = Efieldi(n_levels);
    es_field(:, 1, :) = Pfieldi(:, n_levels);
    es_last_dagger = (Pfieldi(:, n_levels))';
    for fieldi = 2:Nfield
        Vfieldi = V - fieldv(fieldi)*miu;
        Hfieldi = diag(Vfieldi) + Kx;
        [Pfieldi, D] = eig(Hfieldi);
        Efieldi = diag(D);
        es_last_Pfieldi = es_last_dagger*Pfieldi;
        es_last_proj = es_last_Pfieldi.*conj(es_last_Pfieldi);
        [~, max_proj_i] = max(es_last_proj, [], 2);
        Efield(:, fieldi) = Efieldi(max_proj_i);
        es_field(:, fieldi, :) = Pfieldi(:, max_proj_i);
        es_last_dagger = (Pfieldi(:, max_proj_i))';
    end
end