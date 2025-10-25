function U = SchrPWCcheb(Hop, u0, params, Edomain, T, Nt, Ncheb)
% The function performs the propagation of the Schroedinder equation for a
% piecewise constant Hamiltonian. It is assumed that the Hamiltonian
% remains constant within equal time intervals. A Chebyshev algorithm is
% applied for the computation of the function of matrix. Does not apply for
% a non-Hermitian Hamiltonian with complex eigenvalues.
% Input:
% Hop: A function handle of the form @(v, params); represents the operation
% of the Hamiltonian on a vector v. The Hamiltonian depends on a set of
% parametes, specified by a column vector params. Note that one of the
% parametes can be time.
% u0: The initial state vector
% params: A set of time-dependent parameters, which are specified for each 
% time interval. The row index specifies the parametes, and the column
% index specifies the time interval.
% Edomain: The Hamiltonian eigenvalue domain, specified as [minE, maxE].
% T: The final time
% Nt: The number of time-intervals
% Ncheb: The number of Chebyshev points for the computation of the function
% of matrix
% Output:
% U: The state vector evaluated at the initial time and after each time
% interval. The state in different time points is represented by different
% columns.
    dt = T/Nt;
    dim = length(u0);
    U = zeros(dim, Nt + 1);
    U(:, 1) = u0;
    % Computation of the Chebychev coefficients:
    Ccheb = chebc(@(H) exp(-1i*H*dt), Edomain(1), Edomain(2), Ncheb).';
    for ti = 1:Nt
        Hopcheb_ts = @(v) (2*Hop(v, params(:, ti)) - (Edomain(1) + Edomain(2))*v)/(Edomain(2) - Edomain(1));
        v1 = U(:, ti);
        v2 = Hopcheb_ts(U(:, ti));
        U(:, ti + 1) = v1*Ccheb(1) + v2*Ccheb(2);
        for k = 3:Ncheb
            vk = 2*Hopcheb_ts(v2) - v1;
            U(:, ti + 1) = U(:, ti + 1) + vk*Ccheb(k);
            v1 = v2;
            v2 = vk;
        end    
    end
end