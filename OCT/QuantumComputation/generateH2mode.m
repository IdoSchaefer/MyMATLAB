load HTLSs_harmonic
H2modes = multi_kron({I3s, I3s, Manhar, I3s}) + multi_kron({I3s, I3s, I3s, Manhar}) + multi_kron({50*H03e, I3s, I3s, I3s}) +...
    multi_kron({I3s, a3s, adag3s, I3s}) + multi_kron({I3s, adag3s, a3s, I3s})...
    + multi_kron({I3s, a3s, I3s, adag3s}) + multi_kron({I3s, adag3s, I3s, a3s})...
    + multi_kron({a3s, I3s, adag3s, I3s}) + multi_kron({adag3s, I3s, a3s, I3s})...
    + multi_kron({a3s, I3s, I3s, adag3s}) + multi_kron({adag3s, I3s, I3s, a3s});
H2modesu = qubit_excitationsH(H2modes, [3 3 3]);
H2modes1c = multi_kron({I3s, I3s, H03e, I3s});
H2modes1cu = qubit_excitationsH(H2modes1c, [3 3 3]);
H2modes2c = multi_kron({I3s, I3s, I3s, H03e});
H2modes2cu = qubit_excitationsH(H2modes2c, [3 3 3]);
%N2modesu = size(H2modesu, 1);
[Hoperations2modes, fcouplingOp2modes] = Hmats2Hops2(H2modesu, H2modes1cu, H2modes2cu);
%u02modes = zeros(N2modesu, 1);
