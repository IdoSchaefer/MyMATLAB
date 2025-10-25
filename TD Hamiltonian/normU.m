function nU = normU(U)
    nU = sqrt(sum(conj(U).*U));
%     szU = size(U);
%     Nt = szU(2);
%     nU = zeros(1, Nt);
%     for ti = 1:Nt
%         nU(ti) = norm(U(:, ti));
%     end
end