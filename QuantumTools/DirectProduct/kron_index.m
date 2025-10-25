function new_index = kron_index(dims, indices)
% The function computes the index in a direct product vector space from the
% corresponding set of indices in the individual parent spaces. 
% dims: A row vector; represents the dimensions of each of the individual spaces, excluding the first
% component from the left, which is unnecessary for the computation.
% indices: The indices in the individual parent spaces; different columns
% represent the different spaces. It is possible to insert several sets of 
% indices in several different rows.
    Nmats = size(indices, 2);
    dims_prods = cumprod(dims, 'reverse');
    new_index = sum(dims_prods.*(indices(:, 1:(Nmats - 1)) - 1), 2) + indices(:, Nmats);
end