function [Hexcitations, Nsingle, Ndouble, singlei, doublei] = qubit_excitationsH(Hkron, dims)
% The function returns a block diagonal sparse matrix which contains the
% sub-Hamiltonians for a system which is excitation presertive, for the
% initial states in the sub-product-space of two qubits:
% |00>, |01>, |10>, |11>.
% Accordingly, there are 4 blocks, of 0, 1, 1, 2 excitations, respectively. 
% Input:
% Hkron: The Hamiltonian in the direct product space
% dims: A row vector; represents the dimensions of each of the individual spaces, excluding the first
% component from the left, which is unnecessary for the computation.
    singlei = excitation_kroni(1, dims);
    Nsingle = length(singlei);
    doublei = excitation_kroni(2, dims);
    Ndouble = length(doublei);
    Ntotal = 1 + Nsingle + Ndouble;
    Hexcitations = sparse(Ntotal, Ntotal);
    Hexcitations(1, 1) = Hkron(1, 1);
    Hexcitations(2:(Nsingle + 1), 2:(Nsingle + 1)) = Hkron(singlei, singlei);
    Hexcitations((Nsingle + 2):(2*Nsingle + 1), (Nsingle + 2):(2*Nsingle + 1)) = Hkron(singlei, singlei);
    Hexcitations((2*Nsingle + 2):(2*Nsingle + 1 + Ndouble), (2*Nsingle + 2):(2*Nsingle + 1 + Ndouble)) = Hkron(doublei, doublei);
end