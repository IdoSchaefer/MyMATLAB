function result = RK4interp(allx, allt, tresult)
% The function performs a Newton interpolation for a time point, between the
% solution points of ode45.
    Nt = length(allt);
    tgi = ngrater(allt, tresult);
    if allt(tgi) == tresult
        result = allx(:, tgi);
    else    
        if tgi == 2
            spi = 1:4;
        elseif tgi == Nt
            spi = (Nt-3):Nt;
        else
            spi = (tgi-2):(tgi+1);
        end
        result = NewtonIpln4(allt(spi), allx(:, spi), tresult);
    end
end