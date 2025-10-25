function kron_indices = excitation_kroni(Nexcitations, dims)
% The function computes the indices in the direct product space of all
% partition combinations of a given number of excitations into several
% sub-systems.
% Nexcitations: The number of excitations
% dims: A row vector; represents the dimensions of each of the individual spaces, excluding the first
% component from the left, which is unnecessary for the computation.
    Nsystems = length(dims) + 1;
    [~, parent_indices] = excitation_comb(Nexcitations, Nsystems);
    kron_indices = kron_index(dims, parent_indices);
end