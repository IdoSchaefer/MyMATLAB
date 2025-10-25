% Original function
x_i = [0, 1, 3, 4, 7];
f = [1, 3, 49, 129, 813];
% point to evaluate
x = 0.3;

% Setting dimensions: column for x_i and f, row for x
x_i = x_i(:); f = f(:);
x = x(:)';
n = length(x_i);

% The coefficients for the polynomial
% TODO:
%       1. Check the performance for lagre x and f.
%       2. Try with sparse matrices.
DX = tril(bsxfun(@minus,x_i,x_i'));
PD = cumprod([ones(n,1),DX(:,1:end-1)],2);
a = PD\f;
disp('The polynomial coefficients:')
disp(a)

% Evaluating the polynomial at the required points:
% There is an option for multiple points at once
if isscalar(x)
   p_x = a'*[1;cumprod(x-x_i(1:end-1))];
else
   p_x = a'*[ones(1,length(x));cumprod(bsxfun(@minus, x, x_i(1:end-1)) )];
end

disp(' ')
disp('The evaluated polynomial at the desired points:')
disp(p_x)