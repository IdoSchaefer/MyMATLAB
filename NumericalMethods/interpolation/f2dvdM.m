function transM = f2dvdM(z)
% The function computes the transpose of the tranformation matrix from the
% function values to the divided differences.
% The divded differences are given by:
% polcoef = fz*transM,
% where fz represents the function values. Different sampling points are
% represented by different columns.
% The interpolated function can be a vector function.
% Input:
% z: A row vector which contains the sampling points.
    dz1 = bsxfun(@minus, z.', z);
    dz1(dz1==0) = 1;
    prod_dz = cumprod(dz1, 2);
    transM = triu(prod_dz);
    transM(transM~=0) = 1./transM(transM~=0);
end