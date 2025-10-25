function capacity = get_capacity(sp, testp)
% Computation of the capacity of the interpolation domain.
% sp: Vector of the sampling points
% testp: The test point
    capacity = 1;
    sp_comp = sp(sp ~= testp);
    Nsp = length(sp_comp);
    for zi = 1:Nsp
        capacity = capacity*abs(testp - sp_comp(zi))^(1/Nsp);
    end
end
