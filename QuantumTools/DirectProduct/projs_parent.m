function [sq_psi_parent, proj_vecs] = projs_parent(psi, dims, i_parent_space)
    [NHilbert, Ntpoints] = size(psi);
    Nspaces = length(dims);
    sq_psi_parent = zeros(dims(i_parent_space), Ntpoints);
    construction_vecs = cell(1, Nspaces);
    ispaces_ex_parent = 1:Nspaces;
    ispaces_ex_parent(i_parent_space) = [];
    for ispaces = ispaces_ex_parent
        construction_vecs{ispaces} = ones(dims(ispaces), 1);
    end
    construction_vecs{i_parent_space} = sparse(dims(i_parent_space), 1);
    proj_vecs = zeros(NHilbert, dims(i_parent_space));
    for istate = 1:dims(i_parent_space)
        construction_vecs{i_parent_space}(istate) = 1;
        proj_vecs(:, istate) = multi_kron(construction_vecs);
        sq_psi_parent(istate, :) = sum(conj(psi).*proj_vecs(:, istate).*psi);
        construction_vecs{i_parent_space}(istate) = 0;
    end
end