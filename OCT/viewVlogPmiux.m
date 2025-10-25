function viewVlogPmiux(U, Vf, miuxfun, field, x, dt)
    if nargin<6
        dt = 0.1;
    end
    logP = log10(conj(U).*U);
    szU = size(U);
    Nt = szU(2);
    if size(x, 1) == 1
        x = x.';
    end
    field = real(field);
    if length(Vf) == 1
        Vvec = real(Vf(x));
    else
        Vvec = real(Vf);        
        if size(Vvec, 1) == 1
            Vvec = Vvec.';
        end
    end
    if length(miuxfun) == 1
        miux = miuxfun(x);
    else
        miux = miuxfun;
        if size(miux, 1) == 1
            miux = miux.';
        end
    end
    V = (Vvec*ones(1, Nt) - miux*field);
    maxV = max(max(V));
    minV = min(min(V));
    maxP = max(max(logP));
    minP = min(min(logP));
    scalefactor = 5*(maxP - minP)/(maxV - minV);
    V = V*scalefactor;
    figure
    axis([x(1), x(szU(1)), min([minP minV*scalefactor]), max([maxP maxV*scalefactor])])
    line(x,  scalefactor*Vvec, 'color', 'g', 'LineStyle', '--')
    Pcurve = line(x, logP(:, 1));
    Vcurve = line(x, V(:, 1), 'color', 'r');
    for ti = 2:Nt
        pause(dt);
        set(Pcurve, 'Ydata', logP(:, ti))
        set(Vcurve, 'Ydata', V(:, ti))
    end
end