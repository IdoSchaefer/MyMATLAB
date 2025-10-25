function U2q0 = generateU2q0full(dims)
% The function computes the double-qubit vectors, {|00>, |01>, |10>, |11>},
% in the direct product basis, {|mn>}, embedded in a larger direct product
% space. The other subsytems are assumed to be in their ground-state.
% The number of levels in each of the two-qubit subsystems can be larger
% than 2 (e.g. a transmon). Mathematically, It is assumed that the double 
% qubit space is represented by the last two systems from the left in the direct product.
% Input:
% dims: a vector; the dimensions of all the sub-systems in the direct 
% product space
% Output:
% U2q0: A sparse matrix with the double-qubit vectors {|00>, |01>, |10>, |11>}
% represnted by the corresponding columns.
    Ndims = length(dims);
    dim_q2 = dims(Ndims);
    N_Hilbert = prod(dims);
    U2q0 = sparse([1, 2, dim_q2 + 1, dim_q2 + 2], 1:4, ones(1, 4), N_Hilbert, 4);
end