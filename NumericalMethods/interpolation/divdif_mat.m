function polcoef = divdif_mat(z, fz)
% The function computes the divided difference coefficients for a Newton interpolation.
% The program applies also for an interpolated function which returns a vector.
% Input:
% z: A vector which contains the sampling points
% fz: The function values at z. Different sampling points are represented
% by different columns.
% Output:
% polcoef: The coefficients of the Newton basis polynomials for the Newton
% interpolation. The coefficints for the different Newton basis polynomials
% are represented by different columns.
    %[~, Npoints] = size(fz);
    dz1 = bsxfun(@minus, z.', z);
    %dz1 = dz1 + eye(Npoints);
    dz1(dz1==0) = 1;
    prod_dz = cumprod(dz1, 2);
    transM = triu(prod_dz);
    transM(transM~=0) = 1./transM(transM~=0);
    polcoef = fz*transM;
end