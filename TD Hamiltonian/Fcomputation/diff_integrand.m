function result = diff_integrand(m, z, t, n)
    
    gcoef_last = zeros(1, m + 1);
    gcoef_last(m + 1) = 1;
    for deri = 1:n
        gcoef = zeros(1, m + 1);
        for orderi = (m + 1):-1:2
            gcoef(orderi) = gcoef(orderi) - z*gcoef_last(orderi);
            gcoef(orderi - 1) = gcoef(orderi - 1) + (orderi - 1)*gcoef_last(orderi);
        end
        gcoef(1) = gcoef(1) - z*gcoef_last(1);
        gcoef_last = gcoef;
    end
    tpowers = create_tpowerM(t, m);
    result = exp(-z*t).*(gcoef*tpowers);
end

function tpowers = create_tpowerM(t, m)
    Nt = length(t);
    tpowers = zeros(m + 1, Nt);
    tpowers(1, :) = ones(1, Nt);
    for poweri = 1:m
        tpowers(poweri + 1, :) = tpowers(poweri, :).*t;
    end
end