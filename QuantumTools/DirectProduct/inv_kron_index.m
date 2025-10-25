function parent_indices = inv_kron_index(dims, index)
% The inverse function of kron_index; using the index in a direct product
% vector space, the function computes the indices in the individual parent
% spaces.
% dims: A row vector; represents the dimensions of each of the individual spaces, excluding the first
% component from the left, which is unnecessary for the computation.
% index: The index in the direct product space; it is possible to solve the
% problem simultaneously for several indices in the direct product space.
% In this case, index is a column vector.
    Nsets = length(index);
    Nmats = length(dims) + 1;
    parent_indices = zeros(Nsets, Nmats);
    dims_prods = cumprod(dims, 'reverse');
    rem_index = index;
    for dimi = 1:(Nmats - 1)
        % Integer division plus 1:
        parent_indices(:, dimi) = ceil(rem_index./dims_prods(dimi));
        % The remainder:
        rem_index = rem_index - (parent_indices(:, dimi) - 1)*dims_prods(dimi);
    end
    parent_indices(:, Nmats) = rem_index;
end