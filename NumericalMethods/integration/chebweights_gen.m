function W = chebweights_gen(Nsp, domain, result_points)
% The program computes the weights of the function values sampled on
% Chebyshev points that include the boundary, used for integration.
% Exact for a polynomial of degree N-1. The upper limit of integration can
% be any point inside the sampled domain. Weights for several upper limits
% can be computed.
% Input:
% Nsp: The number of Chebyshev sampling points.
% domain: The domain of the sampling points.
% result_points: The upper limits of integration (can be several
% points).
% Output:
% W: The weights of integration for each point in result_points are in the
% corresponding column of W.
    if size(result_points, 1) > 1
        result_points = result_points.';
    end
    lengthD = domain(2) - domain(1);
    Nrp = length(result_points);
    % The result_points mapped to the theta domain, [0 pi]:
    thetap = acos((-2*result_points + domain(1) + domain(2))/lengthD);
    %integT = zeros(Nsp, Nrp);
    n = (0:(Nsp - 1)).';
    % The result of integration over the Chebyshev polynomials, up to a
    % constant factor of 0.5:
    integT = -(cos((n + 1)*thetap) - 1)./(4*(n + 1)*ones(1, Nrp));
    integT([1, 3:Nsp], :) = integT([1, 3:Nsp], :) + ...
        (cos((n([1, 3:Nsp]) - 1)*thetap) - 1)./(4*(n([1, 3:Nsp]) - 1)*ones(1, Nrp));
    W = lengthD*chebcbM(integT);
end