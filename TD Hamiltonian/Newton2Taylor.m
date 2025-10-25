function Ctaylor = Newton2Taylor(x, Cnewton)
% The function computes the coefficients of the form of a Taylor expansion,
% of the interpolation polynom, from the devided differnces.
% x is the vector of the interpolation points.
% Cnewton contains the devided difference vectors in its columns.
    [dim, NC] = size(Cnewton);
    Q = r2Taylor(x);
%    Q = r2Taylor1(x);
    Ctaylor = zeros(dim, NC);
    Ctaylor(:, 1) = Cnewton(:, 1);
    for Newtoni = 2:NC
        for Taylori = 1:Newtoni
            Ctaylor(:, Taylori) = Ctaylor(:, Taylori) + Cnewton(:, Newtoni)*Q(Newtoni, Taylori);
        end
    end
end    