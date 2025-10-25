function v = exact_expAv(A, v0)
    [P, D] = eig(A);
    expD = diag(exp(diag(D)));
    v = P*expD*inv(P)*v0;
end
