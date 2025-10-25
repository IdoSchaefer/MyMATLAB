function [Vopt, optval, meanov] = optabsV(fopt, fguess, NV, xdomain, kdomain, Nk, Nguess, tol, penalty)
% The function finds an optimal absorbing potential Vopt, after repeating
% the optimization procedure for Nguess different guesses for the
% potential.
% Input:
% fopt: A function handle of the form @(Vopt, xdomain, kdomain, Nk, penalty).
% Returns the potential dependent value to be minimized.
% fguess: A function handle of the form @(NV). Returns the guess potential
% column vector. The real part is represented by the first NV/2 terms, and the
% imaginary part by the last NV/2 terms.
% NV: The number of potential elements, including both the real and the imaginary
% elements.
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nk: The number of equally spaced points in the k domain (including the
% boundaries).
% tol: The tolerance of the optimization process
% penalty: Penalty factor for the potential size.
% Output:
% Vopt: The optimal potential found in the process.
% optval: The optimality parameter value of Vopt.
% meanov: The mean of optval for all guesses.
    defaultop = optimset('fminunc');
    options = optimset(defaultop, 'Display', 'off', 'LargeScale', 'off', 'TolFun', tol, 'TolX', tol);
    foptxk = @(V) fopt(V, xdomain, kdomain, Nk, penalty);
%     defaultop = optimset('fmincon');
%     options = optimset(defaultop, 'Display', 'off', 'TolFun', tol, 'TolX', tol, 'Algorithm', 'interior-point');
%     [Vopt, optval] = fmincon(foptxk, fguess(NV), [], [], [], [], -Inf*ones(NV, 1), [ones(NV/2, 1)*Inf; zeros(NV/2, 1)], [], options);
    [Vopt, optval] = fminunc(foptxk, fguess(NV), options);
    meanov = optval;
    for iguess = 2:Nguess
        [Voptnew, optvalnew] = fminunc(foptxk, fguess(NV), options);
        if optvalnew < optval
            Vopt = Voptnew;
            optval = optvalnew;
        end
        meanov = meanov + optvalnew;
    end
    meanov = meanov/Nguess;
end