function polcoef = divdif_morag(z, fz)
% The function computes the divided difference coefficients for a Newton interpolation.
% The routine is based on a divided difference table, where each
% coefficient is given by the last term of a new diagonal in the table.
% The program applies also for an interpolated function which returns a vector.
% Input:
% z: A vector which contains the sampling points
% fz: The function values at z. Different sampling points are represented
% by different rows.
% Output:
% polcoef: The coefficients of the Newton basis polynomials for the Newton
% interpolation. The coefficients for the different Newton basis polynomials
% are represented by different rows.
%     [~, Npoints] = size(fz);
%     diff_z = tril(bsxfun(@minus, z.', z(1:(Npoints - 1))));
%     polyM = cumprod([ones(Npoints, 1), diff_z], 2);
%     polcoef = (polyM\fz.').';
    [Npoints, ~] = size(fz);
    diff_z = tril(bsxfun(@minus, z, z(1:(Npoints - 1)).'));
    polyM = cumprod([ones(Npoints, 1), diff_z], 2);
    polcoef = (polyM\fz);
end