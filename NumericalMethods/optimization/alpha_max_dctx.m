function max_alpha = alpha_max_dctx(x0, direction, max_dctx, dctfactor)
% The function returns the maximal allowed alpha value for the line-search
% procedure line_search1, for a given maximal allowed value of the 
% magnitude of the dctI of the optimization variable vector. It is possible to set a uniform
% value for all the transformed x (the infinity norm of the dctI of x), or to specify the
% value for each of the components of the transformed x.
% x0: The initial point of the search, i.e. the solution from the previous
% iteration
% direction: The direction of search
% max_dctx: The value of the maximal allowed magnitude of the tranformed x for each of the transformed x
% components; if all the components are restricted to a uniform maximal
% magnitude, max_dctx may be a positive scalar. If there are different maximal
% magnitudes for the different components, max_dctx is a vector of the
% dimension of x0.
% dctfactor: The conversion factor between the frequency and time domain.
    dct_direction = dctI(direction)/dctfactor;
    dct_x0 = dctI(x0)/dctfactor;
    max_alpha = min((max_dctx - sign(dct_direction).*dct_x0)./abs(dct_direction));
end