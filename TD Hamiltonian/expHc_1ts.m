function u = expHc_1ts(K, V, u0, leftb, rightb, Ncheb, Ccheb)
% The function computes the next state vector u, from the previous one, u0.
% V is the potential energy. It's a vector in the x domain.
% K is the kinetic energy vector in the p domain.

%    Hpsicheb = @(K, V, psi) (2*Hpsi(K, V, psi) - (leftb + rightb)*psi)/(rightb - leftb);
    v1 = u0;
    v2 = Hpsicheb(K, V, u0, leftb, rightb);
    u = v1*Ccheb(1) + v2*Ccheb(2);
    for k = 3:Ncheb
        vk = 2*Hpsicheb(K, V, v2, leftb, rightb) - v1;
        u = u + vk*Ccheb(k);
        v1 = v2;
        v2 = vk;
    end
end