function v_out = fine2coarse(v_in)
% The function gets a vector which represents a two dimensional space, and
% returns a vector in a coarse grid, with half of the grid points in each
% dimension.
    Nv = sqrt(length(v_in));
    % Note that Nv == N - 1.
    ones_new = ones(floor(Nv/2), 1);
    inject_indices = kron((Nv:2*Nv:Nv*(Nv - 1)).', ones_new) + kron(ones_new, (2:2:Nv).');
    v_out = v_in(inject_indices);
end