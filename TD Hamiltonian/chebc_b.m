function c = chebc_b(f, leftb, rightb, N)
% The function returns the Chebychev coefficients of the function f in a given domain.
% f is a function handle. leftb and rightb are the boundaries of the
% domain. N is the number of Chebychev coefficients.
  %  xcheb = cos(((1:N)*2 - 1)*pi/(2*N));
    xcheb=cos(((1:N) - 1)*pi/(N-1));
    x = 0.5*(xcheb*(rightb - leftb) + rightb + leftb);
%     c = dct(f(x))/sqrt(N);
%     c(2:N) = c(2:N)*sqrt(2);
    f=f(x)';
    c=fft([f; flipud(f(2:N-1))]);    % Extend and compute fft
    c(1)=.5*c(1); c(N)=.5*c(N);
    c=c/(N-1);
    c=c(1:N);
end