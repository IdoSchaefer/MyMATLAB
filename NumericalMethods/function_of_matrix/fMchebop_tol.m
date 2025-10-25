function U = fMchebop_tol(Op, leftb, rightb, f, u0, T, Nt, maxNcheb, tol)
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
    Ccheb = zeros(maxNcheb, Nt);
    for ti = 1:Nt
        ft = @(x) f(x*ti*dt);
        Ccheb(:, ti) = chebc(ft, leftb, rightb, maxNcheb).';
    end
    Opcheb = @(v) (2*Op(v) - (leftb + rightb)*v)/(rightb - leftb);
    v1 = u0;
    v2 = Opcheb(u0);
    U(:, 2:Nt+1) = v1*Ccheb(1,:) + v2*Ccheb(2, :);
    k = 3;
    alltE = ones(1, Nt + 1);
    alltE(1) = 0;
    ti_compute = true(1, Nt + 1);
    ti_compute(1) = false;
    while k<=maxNcheb && max(ti_compute)
        vk = 2*Opcheb(v2) - v1;       
        U(:, ti_compute) = U(:, ti_compute) + vk*Ccheb(k, ti_compute(2:(Nt + 1)));      
        alltE(ti_compute) = norm(vk)*abs(Ccheb(k, ti_compute(2:(Nt + 1))))./(sqrt(sum(U(:, ti_compute).*conj(U(:, ti_compute)))) + 10*eps);
        ti_compute(ti_compute) = alltE(ti_compute)>tol;
        v1 = v2;
        v2 = vk;
        k = k + 1;
    end
    if k>maxNcheb
        display('The program has failrd to achieve the desired tolerance.')
    end
    k
    [maxE, timaxE] = max(alltE)
end