function [H, H1, H2, Hu, H1u, H2u, Hoperations2modes, fcouplingOp2modes] = generateH2mode_gen(alpha, g1, g2, Delta_nu)
    I3s = sparse(eye(3));
    a3s = spdiags(sqrt((0:2).'), 1, 3, 3);
    adag3s = spdiags(sqrt((1:3).'), -1, 3, 3);
    Manhar = spdiags([0; 0; -alpha], 0, 3, 3);
    H3harmonic = spdiags((0:2).', 0, 3, 3);
    H = multi_kron({I3s, I3s, Manhar, I3s}) + multi_kron({I3s, I3s, I3s, Manhar}) + multi_kron({Delta_nu*H3harmonic, I3s, I3s, I3s}) +...
        g2*(multi_kron({I3s, a3s, adag3s, I3s}) + multi_kron({I3s, adag3s, a3s, I3s})...
        + multi_kron({I3s, a3s, I3s, adag3s}) + multi_kron({I3s, adag3s, I3s, a3s}))...
        + g1*(multi_kron({a3s, I3s, adag3s, I3s}) + multi_kron({adag3s, I3s, a3s, I3s})...
        + multi_kron({a3s, I3s, I3s, adag3s}) + multi_kron({adag3s, I3s, I3s, a3s}));
    Hu = qubit_excitationsH(H, [3 3 3]);
    H1 = multi_kron({I3s, I3s, H3harmonic, I3s});
    H1u = qubit_excitationsH(H1, [3 3 3]);
    H2 = multi_kron({I3s, I3s, I3s, H3harmonic});
    H2u = qubit_excitationsH(H2, [3 3 3]);
    %N2modesu = size(H2modesu, 1);
    [Hoperations2modes, fcouplingOp2modes] = Hmats2Hops2(Hu, H1u, H2u);
end
