function U = fMkr(M, f, u0, T, Nt, Nkr, tol)
% The program computes: u(t) = f(M*t)*u0, where M is a matrix, and u0 is a vector. The
% result is computed at Nt equally spaced time points, where T is the maximal
% time. U contains the vector results, in different columns for different
% time points.
% The "Restarted Krylov" algorithm is used.
% f is a function handle.
% Note: The algorithm for computation of the Leja points is unefficient due
% to memory management issues. There exists an enhanced version, Lejap1.m.
% This should be fixed.
    dt = T/Nt;
    dim = length(u0);
%     e1 = zeros(Nkr, 1);
%     e1(1) = 1;
%     I = eye(Nkr);
    U = zeros(dim, Nt + 1);
    U(:, 1) = u0;
    r = u0;
%    r_kr = zeros(Nkr + 1, 1);
%     normr0 = norm(u0);
%     normr = normr0;
    normr = norm(u0);
%    maxNcycles = ceil(dim/Nkr);
    maxNcycles = 100;
    allz = zeros(Nkr, maxNcycles);
    lastdiag = zeros(maxNcycles*Nkr, Nt);
% It might take too much memory. In this case, change maxNcycles.
    if nargin<7
        tol = 1e-5;
%         if nargin<6
%             Nkr = 10;
%         end
    end
    capacity = 0;
    ncycles = 0;
    er = tol + 1;
    while er > tol && ncycles<maxNcycles
        ncycles = ncycles + 1;
        [V, H] = createKr(M, r, Nkr);
%        allz(:, ncycles) = (eig(H(1:Nkr, :)));
%         if ncycles == 1
%             meanz = sum(allz(:, ncycles))/Nkr;
%             newz = allz(allz(:, ncycles) ~= meanz, ncycles);
%             Nnewz = length(newz);
%             for zi = 1:Nnewz
%                 capacity = capacity*abs(meanz - newz(zi))^(1/Nnewz);
%             end
%         end
        [allz(:, ncycles), capacity] = leja_ev(H(1:Nkr, :), allz(:, 1:(ncycles-1)), capacity, Nkr);
%         [P, D] = eig(H(1:Nkr, :));
%         eigval = diag(D);
        [allt_dvd, lastdiag] = newdivdif(allz(:, 1:ncycles), f, lastdiag, capacity, Nkr, ncycles, Nt, dt);
%            U(:, ti + 1) = U(:, ti + 1) + normr*V(:, 1:Nkr)*P*(f(eigval*ti*dt).*(P\e1));
%         for ti = 1:Nt            
%             miu = normr*e1;
%             U(:, ti + 1) = U(:, ti + 1) + allt_dvd(1, ti)*miu;
%             for zi = 1:(Nkr - 1)
%                 miu = (H(1:Nkr, :) - allz(ncycles, zi)*I)*miu;
%                 U(:, ti + 1) = U(:, ti + 1) + allt_dvd(zi + 1, ti)*miu;
%             end
%         end          
%         miu = normr*e1;
%         miu = zeros(Nkr, 1);
%         miu(1) = normr;
%         Ukr = miu*allt_dvd(1, :);
%         for zi = 1:(Nkr - 1)
%             miu = H(1:Nkr, :)*miu - allz(zi, ncycles)*miu;
%             Ukr = Ukr + miu*allt_dvd(zi + 1, :);
%         end
        r_kr = zeros(Nkr + 1, 1);
        r_kr(1) = normr;
        Ukr = zeros(Nkr, Nt);
        for zi = 1:Nkr
            Ukr = Ukr + r_kr(1:Nkr)*allt_dvd(zi, :);
            r_kr(1:(zi + 1)) = (H(1:(zi + 1), 1:zi)*r_kr(1:zi) - allz(zi, ncycles)*r_kr(1:zi+1))/capacity;
% Note all the values of r_kr(>zi) are 0, because of the upper
% triangular form of the Hessenberg matrix H (excluding the first row).
% The Hessenberg matrix and its eigenvalues are divided by the capacity of
% the domain.
        end
        U(:, 2:(Nt + 1)) = U(:, 2:(Nt + 1)) + V(:, 1:Nkr)*Ukr;
%         r_kr(1:Nkr) = normr*(H(1:Nkr, :) - eigval(Nkr)*I)*miu;
%         r_kr(Nkr + 1) = normr*H(Nkr + 1, Nkr)*miu(Nkr);
%         r_kr(1:Nkr) = H(1:Nkr, :)*miu - allz(Nkr, ncycles)*miu;
%         r_kr(Nkr + 1) = H(Nkr + 1, Nkr)*miu(Nkr);
        r = V*r_kr;
        normr = norm(r_kr);
        er = max(abs(allt_dvd(Nkr, :)))*normr;
    end
    ncycles
    if ncycles==maxNcycles && er>tol
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    er
end


%%%% This function can be improved. See fMkrtol.
function [polcoef, diagonal] = newdivdif(allz, f, lastdiag, capacity, Nkr, ncycles, Nt, dt)
    polcoef = zeros(Nkr, Nt);
    Npoints = ncycles*Nkr;
    Noldp = Npoints - Nkr;
%    allzv = allz(1:Npoints);
%    newz = allz(ncycles, :);
    diagonal = lastdiag;
    for ti = 1:Nt
        for coefi = (Noldp + 1):Npoints
            diagonal(coefi, ti) = f(allz(coefi)*ti*dt);
            % Note that we treat allz as a vector, and not a matrix.
            for dtermi = coefi-1:-1:1
                diagonal(dtermi, ti) = (diagonal(dtermi + 1, ti) - diagonal(dtermi, ti))*capacity./(allz(coefi) - allz(dtermi));
% The value of diagonal(:, dtermi) is from the previous diagonal, above diagonal(:, dterm + 1).
% All the sampling points are divided by the capacity.
            end
            polcoef(coefi - Noldp, ti) = diagonal(1, ti);
        end
    end
end

function [ev_out, capacity] = leja_ev(M, old_p, capacity, dimM)
% The function finds all the eigenvalues of the sub matrices of a matrix M.
% Then, it returns the dimM leja points from all the eigenvalues, where dimM
% is the dimension of the matrix.
% old_p: The previous points, for the computation of the new leja points.
% capacity: the capacity of the eigenvalue domain.
    new_ev = zeros(dimM*(dimM + 1)/2, 1);
    % The length of new_ev is the sum of the series: 1+2+3+...+dimM
    evi = 1;
    for dimsub = 1:dimM
        new_ev(evi:(evi + dimsub - 1)) = eig(M(1:dimsub, 1:dimsub));
        evi = evi + dimsub;
    end
    if capacity == 0
        capacity = 1;
        meanz = sum(new_ev)/(evi - 1);
        newz = new_ev(new_ev ~= meanz);
        Nnewz = length(newz);
        for zi = 1:Nnewz
            capacity = capacity*abs(meanz - newz(zi))^(1/Nnewz);
        end
    end
    ev_out = Lejap(old_p(:), new_ev, dimM, capacity);
end