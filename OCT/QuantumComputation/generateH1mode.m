load HTLSs_harmonic
H3harcs10nl = multi_kron({I3s, Manhar, I3s}) + multi_kron({I3s, I3s, Manhar}) + multi_kron({a3s, adag3s, I3s}) + multi_kron({adag3s, a3s, I3s})...
    + multi_kron({a3s, I3s, adag3s}) + multi_kron({adag3s, I3s, a3s});
H3harcus10nl = qubit_excitationsH(H3harcs10nl, [3 3]);
H3LSe_1hs = multi_kron({I3s, H03e, I3s});
H3LSe_1hus = qubit_excitationsH(H3LSe_1hs, [3 3]);
H3LSe_2hs = multi_kron({I3s, I3s, H03e});
H3LSe_2hus = qubit_excitationsH(H3LSe_2hs, [3 3]);
[Hoperations10nl, fcouplingOp3u1ce] = Hmats2Hops2(H3harcus10nl, H3LSe_1hus, H3LSe_2hus);



