function viewVP(U, Vf, field, x, dt)
    if nargin<5
        dt = 0.1;
    end
    P = conj(U).*U;
    szU = size(U);
    Nt = szU(2);
    if size(x, 1) == 1
        x = x.';
    end
    field = real(field);
    V = real(Vf(x)*ones(1, Nt) - x*field);
    maxV = max(max(V));
    minV = min(min(V));
    maxP = max(max(P));
    minP = min(min(P));
    scalefactor = 5*(maxP - minP)/(maxV - minV);
    V = V*scalefactor;
    figure
    axis([x(1), x(szU(1)), min([minP minV*scalefactor]), max([maxP maxV*scalefactor])])
    line(x,  scalefactor*Vf(x), 'color', 'g', 'LineStyle', '--')
    Pcurve = line(x, P(:, 1));
    Vcurve = line(x, V(:, 1), 'color', 'r');
    for ti = 2:Nt
        pause(dt);
        set(Pcurve, 'Ydata', P(:, ti))
        set(Vcurve, 'Ydata', V(:, ti))
    end
end