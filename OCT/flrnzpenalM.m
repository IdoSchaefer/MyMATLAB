function beta = flrnzpenalM(hight, shift, penal, T, dt, Nt_ts, W)

    tcheb = -cos(((1:(Nt_ts - 1)) - 1)*pi/(Nt_ts-1));
    % Note that here, the right boundary is not included - there are
    % Nt_ts-1 points.
    t_ts = 0.5*(tcheb+1)*dt;
    if nargin<7
        W = pi/dt;
    end
    Nt = T/dt;
%     timegrid = zeros(1, 3*Nt*(Nt_ts - 1) + 1);
%     for ti = 1:(Nt_ts - 1)
%         timegrid(ti:(Nt_ts - 1):(3*Nt*(Nt_ts - 1) + 1)) = (-T + t_ts(ti)):dt:2*T;
%     end
    timegrid = zeros(1, Nt*(Nt_ts - 1) + 1);
    for ti = 1:(Nt_ts - 1)
        timegrid(ti:(Nt_ts - 1):(Nt*(Nt_ts - 1) + 1)) = t_ts(ti):dt:T;
    end
    betataufun = @(tau) pi*penal./tau.*(2*hight*shift./tau + sin(W*tau).*(1/hight + ...
        hight*((W - shift)^2*ones(size(tau)) - 2./tau.^2)) + 2*hight./tau.*cos(W*tau)*(W - shift));
    tM = timegrid.'*ones(1, Nt*(Nt_ts - 1) + 1);
    tpM = ones(Nt*(Nt_ts - 1) + 1, 1)*timegrid;
    plustM = tM + tpM;
    minustM = tM - tpM;
    clear tM tpM;
    plusbeta = zeros(Nt*(Nt_ts - 1) + 1);
    plusbeta(plustM == 0) = pi*penal*W*(1/hight + hight*(shift^2 + W^2/3 - W*shift));
    plusbeta(plustM ~= 0) = betataufun(plustM(plustM ~= 0));
    clear plustM;
    minusbeta = zeros(Nt*(Nt_ts - 1) + 1);
    minusbeta(minustM == 0) = pi*penal*W*(1/hight + hight*(shift^2 + W^2/3 - W*shift));
    minusbeta(minustM ~= 0) = betataufun(minustM(minustM ~= 0));
    clear minustM;
    beta = (plusbeta + minusbeta)/2;
end
    