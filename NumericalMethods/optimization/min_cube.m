function minimizer = min_cube(a, b, fa, Dfa, fb, Dfb, interval)
    L = b - a;
    Dz_fa = Dfa*L;
    Dz_fb = Dfb*L;
    Delta = fb - fa;
    c0 = fa;
    c1 = Dz_fa;
    c2 = 3*Delta - 2*Dz_fa - Dz_fb;
    c3 = Dz_fa + Dz_fb - 2*Delta;
    sqrt_expression = c2^2 - 3*c1*c3;
    interval_z = (interval - a)/L;
    if sqrt_expression>0
        zmin = (-c2 + sqrt(sqrt_expression))/(3*c3);
        % According to this specific line-search algorithm, the sign of L
        % is identical to the sign of interval(2)-interval(1), which
        % ensures that interval_z is ordered in increasing values.
        if zmin>interval_z(1) && zmin<interval_z(2)
            % if alpha_min is included in interval:            
            candidates = [interval, a + zmin*L];
            %since alpha_min = a + zmin*L;
            candidates_z = [interval_z, zmin];
        else
            candidates = interval;
            candidates_z = interval_z;
        end
    else
        candidates = interval;
        candidates_z = interval_z;        
    end
    pol_candidates = cube_pol(candidates_z, c0, c1, c2, c3);
    % Actually, this can be made more efficient, since it isn't always necessary
    % to compute the polynomial at the two boundaries. However, it is
    % assumed that the computational time is negligible, and this is
    % much easier to write than a complex system of conditions.
    [~, imin] = min(pol_candidates);
    minimizer = candidates(imin);
end

function result = cube_pol(x, c0, c1, c2, c3)
    result = c0 + c1*x + c2*x.^2 + c3*x.^3;
end
