function gridw = chebgridw(Nt, Nt_ts, dt)
% The function returns a row vector, containig the weights of integration, for
% the points in a grid devided into equal intervals with internal Chebyshev
% structure.
    % Computation of the weights of the Chebyshev points in integration:
    weightcheb = chebweights(Nt_ts, dt);
    % The weights of the whole grid. The weight of the points on the
    % boundary between 2 time steps is multiplied by 2:
    gridw = [weightcheb(1:(Nt_ts - 1)), kron(ones(1, Nt - 1), [weightcheb(1)*2, weightcheb(2:(Nt_ts - 1))]), weightcheb(Nt_ts)];
end