function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn(psi0, H0, Edomain, miu,...
    fguess, filterE, filtermiu, options, T, dt, Nt_ts, Ncheb, tol, maxNiter)

    Nt = T/dt;
    dw = pi/T;
    integw = [dw/2; ones(Nt - 1, 1)*dw; dw/2];
    dctfactor = T/(sqrt(Nt*pi));
    dim = length(psi0);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    tolprop = 1e-3*tol;
    summniterc = 0;
    allfield = zeros(allt_lasti, 1);
    allevmiut = zeros(allt_lasti, 1);
    evmiufil = zeros(allt_lasti + 1, 1);
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
%    Nwnz = length(iEnz);
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
    vfiltermiu = filtermiu(w);
    allpsi = zeros(dim, allt_lasti);
    allmiupsi = zeros(dim, allt_lasti);
    chiT = zeros(dim, 1);
    chiihterm = zeros(dim, allt_lasti);
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
        [allpsi, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield.');
        summniterc = summniterc + mniterc;
        allmiupsi = miu*allpsi;
        for allti = 1:allt_lasti
% The values of the expectation value are supposed to be real:
            allevmiut(allti) = real(allpsi(:, allti)'*allmiupsi(:, allti));
        end
        evmiuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        get_chiihterm(allmiupsi);
        [allchi, ~, mniterc] = solveOCMLSih1(@ihalltchiMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1), chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        J1fun = 0.5*evmiuw.^2.*vfiltermiu.*integw;
        J2fun = -(fieldwnz.^2./vfilterEnz).*integwnz;
        minusJ = -(sum(J1fun) + sum(J2fun));
        for allti = 1:allt_lasti
            chimiupsi(allti) = imag(allchi(:, allt_lasti - allti + 1)'*(miu*allpsi(:, allti)));
        end
        dct_chimiupsi = dctfactor*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights);
        minusgrad = 2*integwnz.*(fieldwnz./vfilterEnz + dct_chimiupsi(iEnz));
    end


    function get_chiihterm(allmiupsi)
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        chiihterm = -allmiupsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1), 0, allt_lasti, allt_lasti);     
    end

end