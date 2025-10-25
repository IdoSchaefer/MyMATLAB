function viewP(U, x, dt, labelx, labely)
% The function views the evolution of the probability of the time-dependent
% state U.
% Input:
% U: Matrix; The time-dependent state evaluated at equally spaced time-points;
% each time-point is represented by a different column.
% x: The grid of U
% dt: The time spacing between the views of consecutive columns of U; the
% default value is 0.1.
% labelx: String; the x label of the plot, where latex is the interpreter;
% the default value is 'x [a.u.]'.
% labely: String; the y label of the plot, where latex is the interpreter;
% the default value is '$\left|\psi(x)\right|^2$'.
    if nargin<3
        dt = 0.1;
    end
    if nargin<4
        labelx = 'x [a.u.]';
    end
    if nargin<5
        labely = '$\left|\psi(x)\right|^2$';
    end
    P = conj(U).*U;
    szU = size(U);
    Nt = szU(2);
    figure
    axis([x(1), x(szU(1)), min(min(P)), max(max(P))])
    xlabel(labelx, 'interpreter', 'latex')
    ylabel(labely, 'interpreter', 'latex')
    Pcurve = line(x, P(:, 1));
    for ti = 2:Nt
        pause(dt);
        set(Pcurve, 'Ydata', P(:, ti))
    end
end