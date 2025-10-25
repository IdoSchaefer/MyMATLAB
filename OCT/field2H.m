function H = field2H(field, E)
    szf = size(field);
    dim = 2;
%    dim = 2*szf(1);
    Nt = szf(2);
    H = zeros(dim, dim, Nt);
    H(1, 1, :) = E(1);
    H(2, 2, :) = E(2);
    for ti = 1:Nt
        H(1, 2) = -conj(field(ti));
        H(2, 1) = -field(ti);
    end
end