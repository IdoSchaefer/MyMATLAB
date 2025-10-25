function ev_out = leja_ev(M, old_p, capacity, dimM)
% The function finds all the eigenvalues of the sub matrices of a matrix M.
% Then, it returns the dimM leja points from all the eigenvalues, when dimM
% is the dimention of the matrix.
% old_p: The previous points, for the computation of the new leja points.
% capacity: the capacity of the eigenvalue domain.
    new_ev = zeros(dimM*(dimM + 1)/2, 1);
    % The length of new_ev is the sum of the series: 1+2+3+...+dimM
    evi = 1;
    for dimsub = 1:dimM
        new_ev(evi:(evi + dimsub - 1)) = eig(M(1:dimsub, 1:dimsub));
        evi = evi + dimsub;
    end
    ev_out = Lejap(old_p(:), new_ev, dimM, capacity);
end