function U = fMchebop(Op, leftb, rightb, f, u0, T, Nt, Ncheb)
% The program computes: u(t) = f(Op*t)*u0, where Op is an operator, and u0 is a vector. The
% result is computed at Nt equally spaced time points, where T is the maximal
% time. U contains the vector results, in different columns for different
% time points.
% The program uses an expansion of Chebychev polynomials for f(Op*t).
% f and Op are function handles.
% leftb and rightb are the (real) boundaries of the domain
% of the eigenvalues of M. Ncheb is the number of terms used to expand
% the function in the Chebychev series.
    dt = T/Nt;
    dim = length(u0);
    U = zeros(dim, Nt + 1);
    U(1:dim, 1) = u0;
    % Computation of Chebychev coefficients in all the time points:
    Ccheb = zeros(Ncheb, Nt);
    for ti = 1:Nt
        ft = @(x) f(x*ti*dt);
        Ccheb(:, ti) = chebc(ft, leftb, rightb, Ncheb).';
    end
    Opcheb = @(v) (2*Op(v) - (leftb + rightb)*v)/(rightb - leftb);
    v1 = u0;
    v2 = Opcheb(u0);
    U(:, 2:Nt+1) = v1*Ccheb(1,:) + v2*Ccheb(2, :);
    for k = 3:Ncheb
        vk = 2*Opcheb(v2) - v1;
        U(:, 2:Nt + 1) = U(:, 2:Nt + 1) + vk*Ccheb(k, :);
        v1 = v2;
        v2 = vk;
    end
end