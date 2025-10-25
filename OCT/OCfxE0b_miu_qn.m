function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(psi0, Vf, m,...
    Edomain, xdomain, miux, fguess, filterE, filtermiu, options, T, dt, Nt_ts, Ncheb, tol, maxNiter)

    Nt = T/dt;
    dw = pi/T;
    integw = [dw/2; ones(Nt - 1, 1)*dw; dw/2];
    dctfactor = T/(sqrt(Nt*pi));
    Nx = length(psi0);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/(2*m);
    miuxM = spdiags(miux, 0, Nx, Nx);
    tolprop = 1e-3*tol;
    summniterc = 0;
    allfield = zeros(allt_lasti, 1);
    allevmiut = zeros(allt_lasti, 1);
    evmiufil = zeros(allt_lasti, 1);
    evmiuw = zeros(1, Nt + 1);
    J1fun = zeros(Nt + 1, 1);
    chimiupsi = zeros(allt_lasti, 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
%     allt_grid = [kron((0:dt:(T - dt)).', ones(Nt_ts - 1, 1)) + kron(ones(Nt, 1), t_ts(1:(Nt_ts - 1))); T];
    w = (0:pi/T:pi/dt).';
    vfilterE = filterE(w);
%    [iEnz, ~, vfilterEnz] = find(vfilterE);
%    iEnz = find(vfilterE>=tolprop*max(vfilterE)/10);
    iEnz = find(vfilterE>=eps*max(vfilterE));
    vfilterEnz = vfilterE(iEnz);
    igweights = chebweights(Nt_ts, 1);
    igweights = [2*igweights(1), igweights(2:(Nt_ts - 1))];
    integwnz = integw(iEnz);
    dctfilterE0 = sqrt(2/pi)*sum(vfilterEnz.*integwnz);
    coswT = ones(Nt + 1, 1);
    coswT(2:2:(Nt + 1)) = -1;
    dctfilterET = sqrt(2/pi)*sum(vfilterEnz.*coswT(iEnz).*integwnz);
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    fieldw = zeros(Nt + 1, 1);
    if length(fguess) == 1
        % if fguess is a function handle:
        for wnzi = iEnz
            fieldw(wnzi) = fguess(w(wnzi));
        end
    else
        % if fguess is a vector:
        fieldw(iEnz) = fguess(iEnz);
    end
    fieldwnz = fieldw(iEnz);
    if sqrt(2/pi)*abs(sum(fieldwnz.*integwnz)) > tolprop*1e-2 || sqrt(2/pi)*abs(sum(fieldwnz.*coswT(iEnz).*integwnz)) > tolprop*1e-2
        fieldwnz = getfieldwcon(fieldwnz);
    end
    vfiltermiu = filtermiu(w);
    allpsi = zeros(Nx, allt_lasti);
    allmiupsi = zeros(Nx, allt_lasti);
    chiT = zeros(Nx, 1);
    chiihterm = zeros(Nx, allt_lasti);
    nprop = 0;
    if isempty(options)
        options = optionsOCqn(tol, maxNiter);
    end
    if isempty(options.invHess0)
        options.invHess0 = diag(vfilterEnz./(2*integwnz));
    end
    [fieldwnz, ~, minusgrad, niter, ~, ~, dif_fieldw, minus_conv, alpha, invHess] = quasiNewton(@Jeval, fieldwnz, options);
    fieldw = fieldw.';
    maxgrad = max(abs(minusgrad));
    conv = -minus_conv;
    relE = norm(dif_fieldw)/norm(fieldwnz);
    niter
    nprop
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = allfield(1:(Nt_ts - 1):allt_lasti).';
    evmiut = allevmiut(1:(Nt_ts - 1):allt_lasti).';
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    beep
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad] = Jeval(fieldwnz)
        fieldw(iEnz) = fieldwnz;
        allfield = dctIintgrid(fieldw, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, mniterc] = solveOCih2(@ihfieldmiux, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield.', miux);
        summniterc = summniterc + mniterc;
        allmiupsi = miuxM*allpsi;
        for allti = 1:allt_lasti
% The values of the expectation value are supposed to be real:
            allevmiut(allti) = real(allpsi(:, allti)'*allmiupsi(:, allti));
        end
        evmiuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        get_chiihterm(allmiupsi);
        [allchi, mniterc] = solveOCih2(@ihalltchimiux, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1).', miux, chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        J1fun = 0.5*evmiuw.^2.*vfiltermiu.*integw;
        J2fun = -(fieldwnz.^2./vfilterEnz).*integwnz;
        minusJ = -(sum(J1fun) + sum(J2fun));
        for allti = 1:allt_lasti
            chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*(miux.*allpsi(:, allti)));
        end
        dct_chimiupsi = dctfactor*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights);
        newfieldw_unc_nz = dct_chimiupsi(iEnz).*vfilterEnz;
        [lambda0, lambdaT] = fieldw_unc2lambda(newfieldw_unc_nz);
        minusgrad = 2*integwnz.*(fieldwnz./vfilterEnz - dct_chimiupsi(iEnz) + lambda0 + lambdaT*coswT(iEnz));
    end

    function [lambda0, lambdaT] = fieldw_unc2lambda(fieldw_unc_nz)
        fieldt_unc0 = sqrt(2/pi)*sum(fieldw_unc_nz.*integwnz);
        fieldt_uncT = sqrt(2/pi)*sum(fieldw_unc_nz.*coswT(iEnz).*integwnz);
        lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
        lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;        
    end

    function fieldw_con = getfieldwcon(fieldw_unc_nz)
    % The function computes a constrained field spectrum, fieldw_con, with 0 boundaries in the
    % time domain, from the unconstrained field spectrum fieldw_unc.
        [lambda0, lambdaT] = fieldw_unc2lambda(fieldw_unc_nz);
        fieldw_con = fieldw_unc_nz - vfilterEnz.*(lambda0 + lambdaT*coswT(iEnz));
    end

    function get_chiihterm(allmiupsi)
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        chiihterm = -allmiupsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1), 0, allt_lasti, allt_lasti);     
    end

end