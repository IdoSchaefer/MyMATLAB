function Vt = field2Vt(field)
    szf = size(field);
    dim = 2;
%    dim = 2*szf(1);
    Nt = szf(2);
    Vt = zeros(dim, dim, Nt);
    for ti = 1:Nt
        Vt(1, 2, ti) = -conj(field(ti));
        Vt(2, 1, ti) = -field(ti);
    end
end