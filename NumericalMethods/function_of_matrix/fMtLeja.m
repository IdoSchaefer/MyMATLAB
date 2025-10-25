function [U, Np_final, maxE]  = fMtLeja(Moperation, evdomain, f, u0, T, Nt, Npcycle, maxNp, tol)
% The program computes: u(t) = f(M*t)*u0, where M is a matrix, and u0 is a vector. The
% result is computed at Nt equally spaced time points, where T is the maximal
% time. U contains the vector results, in different columns for different
% time points.
% The function of matrix is computed by a Newton interpolation in Leja
% points in a length 4 domain.
% The program is suitable for M with real eigenvalues.
% Input:
% Moperation: A function handle that defines the operation of M on a
% vector.
% evdomain: The estimated eigenvalue domain of M.
% f is a function handle.
% T: The final time
% Nt: The number of equally spaced time points (excluding t=0).
% Npcycle: The number of additional sampling points used in each step.
% maxNp: The maximal allowed number of sampling points. Determines also the
% size of the reservoir from which the Leja points are chosen.
% tol: The tolerance of the solution.
% Output:
% Np_final: The final number of sampling points used in the interpolation.
% maxE: The maximal estimated error.
    dt = T/Nt;
    dim = length(u0);
    U = zeros(dim, Nt + 1);
    U(:, 1) = u0;
    % Qu0 = Q(M*t)*u0, where Q is the Newton expansion polynomial.
    Qu0 = u0;
    Dsize = evdomain(2) - evdomain(1);
    D4factor = 4/Dsize;
    evdomain4 = evdomain*D4factor;
% The Leja points are chosen from the variable reservoir:
    reservoir = (evdomain4(1):(4/(maxNp*100)):evdomain4(2)).';
    % An alternative reservoir:
    % reservoir = (rand(maxNp*100 + 1, 1) - 0.5)*4 + (evdomain4(1) + evdomain4(2))/2;
    % is_unused is a logical that tells which points from the reservoir are
    % still unused.
    is_unused = true(maxNp*100 + 1, 1); 
    % The multiplicity of the differences from the previous Leja points:
    difmul = ones(maxNp*100 + 1, 1);
    oldp = zeros(maxNp, 1);
% The new divided differences in all time points:
    allt_dvd = zeros(Npcycle, Nt);
    maxNcycles = floor(maxNp/Npcycle);
    diagonal = zeros(maxNp, Nt);
% It might take too much memory. In such a case, change maxNp, or write
% another program...
% Default tolerance:
    if nargin<9
        tol = 1e-5;
    end
    ncycles = 0;
    alltE = tol + ones(1, Nt);
    not_converged = alltE>tol;
    % not_converged is a logical variable, which tells in which time points to
    % continue the computation, because the desired accuracy hasn't been
    % achieved yet.
    while max(not_converged) && ncycles<maxNcycles
        % max(not_converged) == true if there are still unconverged time
        % points.
        % Getting the new Leja points:
        newp = Lejap;
        ncycles = ncycles + 1;
        % Getting the new coefficients (divided differences) for the Newton
        % interpolation:
        newdivdif;
        % Computing the new terms in the Newton expansion:
        for pi = 1:Npcycle
            U(:, [false, not_converged]) = U(:, [false, not_converged]) + Qu0*allt_dvd(pi, not_converged);
            Qu0 = Moperation(Qu0)*D4factor - newp(pi)*Qu0;
        end
        % This is a measure of the error:
        alltE(not_converged) = abs(allt_dvd(Npcycle, not_converged))*norm(Qu0);
        not_converged(not_converged) = alltE(not_converged)>tol;
    end
    Np_final = ncycles*Npcycle
    if ncycles==maxNcycles && max(not_converged)
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    [maxE, timaxE] = max(alltE)
    
    %%%%%%%% Nested functions: %%%%%%%%
    
    function newp = Lejap
        newi = 1;
        nold = ncycles*Npcycle;
% The case of an empty oldp requires a special treatment for the 1'st
% point - we have to choose it smartly:
        if ncycles == 0        
            oldp(1) = reservoir(1);
%             reservoir = reservoir(2:end);
%             difmul = difmul(2:end);
            is_unused(1) = false;
            difmul(1) = 0;
            newi = newi + 1;
        end
        for newi = newi:Npcycle
            %difmul = difmul.*abs(reservoir - oldp(nold + newi - 1));
            difmul(is_unused) = difmul(is_unused).*abs(reservoir(is_unused) - oldp(nold + newi - 1));
            [stam, maxi] = max(difmul);
            oldp(nold + newi) = reservoir(maxi);
%             reservoir = reservoir([1:(maxi-1), (maxi+1):end]);
%             difmul = difmul([1:(maxi-1), (maxi+1):end]);
            % Using the profiler, I found that the 2 last commented lines took the
            % most of the computational effort in this program. Hence, I
            % use now another, less direct technique, with the logical is_unused.
            difmul(maxi) = 0;
            is_unused(maxi) = false;
        end
        newp = oldp((nold + 1):(nold + Npcycle));
    end

    function newdivdif
        Npoints = ncycles*Npcycle;
        Nplast = Npoints - Npcycle;
        ti_compute = 1:Nt;
        ti_compute = ti_compute(not_converged);
        for coefi = (Nplast + 1):Npoints
            diagonal(coefi, ti_compute) = f(oldp(coefi)/D4factor*ti_compute*dt);
            for dtermi = coefi-1:-1:1
                diagonal(dtermi, ti_compute) = ...
                    (diagonal(dtermi + 1, ti_compute) - diagonal(dtermi, ti_compute))/(oldp(coefi) - oldp(dtermi));
% The value of diagonal(dtermi, ti_compute) is from the previous diagonal, above diagonal(dtermi, ti_compute) in 
% the divided differences table.
            end
            allt_dvd(coefi - Nplast, ti_compute) = diagonal(1, ti_compute);
        end
    end

end