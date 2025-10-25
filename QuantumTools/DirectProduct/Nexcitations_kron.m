function Nex = Nexcitations_kron(dims, indices)
% The function computes the number of excitations in a direct product
% space from the indices of the direct product space. It is assumed that
% the parent system spaces are non-degenerate.
% dims: A row vector; represents the dimensions of each of the individual 
% spaces, excluding the first component from the left, which is unnecessary
% for the computation.
% indices: A column vector; The indices in the direct product space.
    Nex = sum(inv_kron_index(dims, indices) - 1, 2);
end