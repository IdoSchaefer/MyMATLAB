function c = chebcv(fv)
% The function computes the Chebychev coefficients of a
% function sampled at the regular Chebychev points 
% (the ones do not include the boundaries of the domain).
% fv is the vector of the sampled values.
    N = length(fv);
    c = dct(fv)/sqrt(N);
    c(2:N) = c(2:N)*sqrt(2);
end