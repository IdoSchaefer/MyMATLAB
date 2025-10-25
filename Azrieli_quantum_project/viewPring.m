function viewPring(U, phi, r, dt)
% The function views the evolution of the probability of the time-dependent
% state U of a particle on a ring. The probability is plotted on a ring of
% radius r.
% Input:
% U: Matrix; The time-dependent state evaluated at equally spaced time-points;
% each time-point is represented by a different column.
% phi: The grid of U
% r: The radius of the ring, for visualization purposes.
% dt: The time spacing between the views of consecutive columns of U; the
% default value is 0.1.
    if nargin<3
        dt = 0.1;
    end
    phi_ext = [phi; 2*pi];
    P = [conj(U).*U; conj(U(1, :)).*U(1, :)];
    Pplot = r*(1 + P/(max(max(P))*3));
    [Nphi, Nt] = size(U);
    figure
    ax = polaraxes;
    ax.RLim = [0, 4*r/3];
    line(phi_ext, r*ones(Nphi + 1), 'color', 'k', 'linewidth', 3)
    Pcurve = line(phi_ext, Pplot(:, 1), 'color', 'm');
    for ti = 2:Nt
        pause(dt);
        set(Pcurve, 'Ydata', Pplot(:, ti))
    end
end