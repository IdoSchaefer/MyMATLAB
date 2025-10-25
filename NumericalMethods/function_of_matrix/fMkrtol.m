function [U, alltNkr, alltE] = fMkrtol(M, f, u0, T, Nt, tol)
% The program computes: u(t) = f(M*t)*u0, where M is a matrix, and u0 is a vector. The
% result is computed at Nt equally spaced time points, where T is the maximal
% time. U contains the vector results, in different columns for different
% time points.
% The simple Krylov algorithm is used. The Krylov space is enlarged in each
% cycle by one dimension.
% f is a function handle.
    dt = T/Nt;
    dim = length(u0);
    H = zeros(dim, dim - 1);
    V = zeros(dim, dim);
    allt_dvd = zeros(dim, Nt);
    diagonal = zeros(dim, Nt);
    U = zeros(dim, Nt + 1);
    Ukr = zeros(dim - 1, Nt);
    eigval = zeros(dim - 1, 1);
    r_kr = zeros(dim, 1);
    is_not_converged = true(1, Nt);
    condition = true(1, Nt);
    normU = zeros(1, Nt);
    alltE = zeros(1, Nt);
    alltNkr = zeros(1, Nt);
    t = dt:dt:T;
    U(:, 1) = u0;
    nu0 = norm(u0);
    V(:, 1) = u0/nu0;
    vj = 1;
    enlargeKr
    vj = 2;
    while max(is_not_converged) && vj < dim
        enlargeKr
%         if vj == 40
%             keyboard
%         end
        eigval(1:vj) = eig(H(1:vj, 1:vj));
        avgp = sum(eigval(1:vj))/vj;
        samplingp = [eigval(1:vj); avgp];
        capacity = get_capacity;
        sp_capacity1 = samplingp/capacity;
%        get_allt_dvd(samplingp, f(samplingp*t(is_not_converged)));
        get_allt_dvd(sp_capacity1, f(samplingp*t(is_not_converged)));
        %allt_dvd(1:vj, is_not_converged) = devdif(eigval(1:vj), f(t(is_not_converged).'*eigval(1:vj).')).';
        r_kr(1) = nu0;
        r_kr(2:vj) = 0;
        Ukr(1:(vj - 1), is_not_converged) = 0;
        for evi = 1:vj
            Ukr(1:vj, is_not_converged) = Ukr(1:vj, is_not_converged) + r_kr(1:vj)*allt_dvd(evi ,is_not_converged);
%            r_kr(1:(evi + 1)) = H(1:(evi + 1), 1:evi)*r_kr(1:evi) - eigval(evi)*r_kr(1:(evi + 1));
            r_kr(1:(evi + 1)) = H(1:(evi + 1), 1:evi)*r_kr(1:evi)/capacity - sp_capacity1(evi)*r_kr(1:(evi + 1));
% Note all the values of r_kr(>evi) are 0, because of the upper
% triangular form of the Hessenberg matrix H (excluding the first row).
% The Hessenberg matrix and its eigenvalues are divided by the capacity of
% the domain.
        end
        normU(is_not_converged) = sqrt(sum(Ukr(:, is_not_converged).*conj(Ukr(:, is_not_converged))));
        alltE(is_not_converged) = abs(allt_dvd(vj + 1, is_not_converged))*norm(r_kr(1:(vj + 1)))./normU(is_not_converged);
        condition(is_not_converged) = alltE(is_not_converged) > tol;
        alltNkr(is_not_converged & ~condition) = vj;
        is_not_converged(is_not_converged) = condition(is_not_converged);
        vj = vj + 1;
    end
    if vj == dim
        display('It was required to diagonalize the entire matrix for the requested tolerance.')
        [P, D] = eig(M);
        U(:, [false, is_not_converged]) = P*(spdiags(P\u0, 0, dim, dim)*f(diag(D)*t(is_not_converged)));
    end
    U(:, [false, ~is_not_converged]) = V(:, 1:(vj - 1))*Ukr(1:(vj - 1), ~is_not_converged);
    
    %%%%%%%%% Nested functions: %%%%%%%%%%
    function get_allt_dvd(x, fx)
        Npoints = vj + 1;
        allt_dvd(1, is_not_converged) = fx(1, :);
        diagonal(1, is_not_converged) = fx(1, :);
        for coefi = 2:Npoints
            diagonal(coefi, is_not_converged) = fx(coefi, :);
            for dtermi = coefi-1:-1:1
                diagonal(dtermi, is_not_converged) =...
                    (diagonal(dtermi + 1, is_not_converged) - diagonal(dtermi, is_not_converged))./(x(coefi) - x(dtermi));
    % The value of diagonal(dtermi, is_not_converged) is from the previous diagonal, above diagonal(dterm + 1, is_not_converged).
            end
            allt_dvd(coefi, is_not_converged) = diagonal(1, is_not_converged);
        end
    end

    function enlargeKr
        V(:, vj+1) = M*V(:, vj);
        for vi = 1:vj
            H(vi, vj) = V(:, vi)'*V(:, vj+1);
            V(:, vj+1) = V(:, vj+1) - H(vi, vj)*V(:, vi);
        end
        H(vj+1, vj) = norm(V(:, vj+1));
        V(:, vj+1) = V(:, vj+1)/H(vj+1, vj);
    end

    function capacity = get_capacity
        capacity = 1;
        sp_comp = samplingp(samplingp ~= samplingp(vj + 1));
        Nsp = length(sp_comp);
        for zi = 1:Nsp
            capacity = capacity*abs(samplingp(vj + 1) - sp_comp(zi))^(1/Nsp);
        end
    end

end

