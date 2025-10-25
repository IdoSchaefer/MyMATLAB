function Mresult = multi_kron(Mats)
% The function computes the Kronecker product of several matrices.
% Mats: A cell array which contains the matrices to be multiplied.
    NM = length(Mats);
    if NM == 2
        Mresult = kron(Mats{1}, Mats{2});
    elseif NM>2
        Mresult = kron(Mats{1}, multi_kron(Mats(2:NM)));
    else
        Mresult = [];
    end
end
