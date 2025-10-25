function [field, psi, relE, conv, niter] = OCmyRK(psi0, target, Vf, max_x, fguess, Epenal, T, dt, tol)
% fguess is a function handle of the form: @(t).
% Vf is the time independent potential energy. It's a function handle of
% the form: f(x), when x is the row vector of the x grid.
% The space domain is: [-max_x, max_x).
    Nt = T/dt;
    Nx = length(psi0);
    dx = 2*max_x/Nx;
    x = (-max_x:dx:(max_x - dx)).';
    p = (0:(pi/max_x):(2*pi/dx - pi/max_x)).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
%     chi = zeros(Nx, Nt+1);
%     psi = zeros(Nx, Nt+1);
    output = zeros(Nx + 1, Nt + 1);
tic
%     psisol = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
%    [t, psi] = ode45(@OCguess, 0:dt:T, psi0, [], K, Vf, x, fguess);
   [tpsi, psisol] = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
%     psisol = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
%     tpsi = psisol.x;    
   psisol = psisol.';
toc
% an estimation of dt for computation, for all the propagations (of course,
% not nessesarily right):
    dtcomp = min(tpsi(2:end) - tpsi(1:end-1));
%    dtcomp = 4e-4;
    Ntcomp = floor(T/dtcomp) + 1;
% for the time interval to contain an (approximately) integer number of time steps:
    dtcomp = T/Ntcomp;
    psi = zeros(Nx, Ntcomp+1);
    chi = zeros(Nx, Ntcomp+1);
    psi(:, 1) = psi0;
    tic
    for ti = 2:(Ntcomp+1)
        psi(:, ti) = RK4interp(psisol, tpsi, (ti - 1)*dtcomp);
    end
%    psi = deval(psisol, 0:dtcomp:T);
    toc
    Ner = max([floor(Ntcomp/100), 1]);
    field = zeros(1, Ntcomp+1);
    lastfield = zeros(1, Ntcomp+1);
%    psiT = deval(psisol, T);
    psiT = psi(:, end);
    for ti = 1:(Ntcomp+1)
        lastfield(ti) = fguess((ti-1)*dtcomp);
    end
    if nargin>3
        field = lastfield;
    end
%    lastfield = field;
    relE = tol + 1;
%    conv = zeros(2, 1001);
    conv = zeros(1, 1001);
    maxNiter = 100;
    niter = 0;
    while relE>tol && niter<maxNiter
        tic
        overlap = target'*psiT;
        chiT = overlap*target;
        if nargin>3
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dtcomp;
        end
%        chisol = ode45(@OCchi, [T 0], chiT, [], K, Vf, x, psisol, Epenal);
%        [t, chi] = ode45(@OCchi, T:-dt:0, chiT, [], K, Vf, x, psi, Epenal, T, dt);
        chi(:, (Ntcomp+1):-1:1) = RK4OC(@OCchim, [T 0], chiT, -dtcomp, K, Vf, x, Nx, psi(:, (Ntcomp+1):-1:1), Epenal);
%        chi = chi(:, end:-1:1);
%        chi = deval(chisol, 0:dt:T);
        [psi field] = RK4OC(@OCpsim, [0 T], psi0, dtcomp, K, Vf, x, Nx, chi, Epenal);
        psiT = psi(:, Nt + 1);
%        psi = deval(psisol, 0:dt:T);        
%         for ti = 1:(Nt+1)
%             chit = RK4interp(chisol, tchi, (ti - 1)*dt);
%             psit = RK4interp(psisol, tpsi, (ti - 1)*dt);
%             field(ti) = -imag(chit'*(x.*psit))/Epenal;
%            field(ti) = RK4interp(psisol.y(Nx + 1, :), psisol.x, (ti - 1)*dt);
%         end
%         psiT = deval(psisol, T);
%         psiT = psiT(1:Nx);
%         for ti = 1:(Ntcomp+1)
%             field(ti) = -imag(chi(:, ti)'*(x.*psi(:, ti)))/Epenal;
%         end
        field(Ntcomp + 1) = -imag(chi(:, Nt + 1)'*(x.*psi(:, Nt + 1)))/Epenal;
        relE = norm(field(1:Ner:(Ntcomp + 1)) - lastfield(1:Ner:(Ntcomp + 1)))/norm(field(1:Ner:(Ntcomp + 1)));
        lastfield = field;
        niter = niter + 1;
        toc
    end
    for ti = 1:(Nt+1)
        output(:, ti) = RK4interp([psi; field], 0:dtcomp:T, (ti - 1)*dt);
    end
    psi = output(1:Nx, :);
    field = output(Nx + 1, :);
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
%    psi = deval(psisol, 0:dt:T);
    if nargin>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = target'*psiT;
        conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dtcomp;
        conv = conv(:, 1:(niter + 1));
    end
    niter
end

function [U field] = RK4OC(fderiv, t0tf, u0, dt, K, Vf, x, Nx, conjU, Epenal)
% fderiv:  a function handle. The value of the function is the derivative of u with respect to t.
% t0tf:  a vector that contains initial and final time: [t0 tf]
% u0:  a vector of the initial values of the result vector u.
% dt: The time step.
% U:  The result u vectors at all the times.
    t0 = t0tf(1);
    tf = t0tf(2);
    Nt = floor(abs((tf - t0)/dt));
%    dim = length(u0);
    U = zeros(Nx, Nt + 1);
    field = zeros(1, Nt + 1);
    U(:, 1) = u0;
    tiIpln = 1:4;
    conjUhalf = NewtonIpln4(t0 + (tiIpln-1)*dt, conjU(:, tiIpln), t0 + 0.5*dt);
    [U(:, 2) field(1)] = unew(fderiv, u0, K, Vf, x, conjU(:, 1), conjUhalf, conjU(:, 2), Epenal, dt);
    t = t0 + dt;
    for ti = 2:Nt-1
        tiIpln = (ti - 1):(ti + 2);
        conjUhalf = NewtonIpln4(t0 + (tiIpln-1)*dt, conjU(:, tiIpln), t + 0.5*dt);
        [U(:, ti + 1) field(ti)] = unew(fderiv, U(:, ti), K, Vf, x, conjU(:, 1), conjUhalf, conjU(:, 2), Epenal, dt);
 % The k's are estimations of du. The result is composed of a weighted
 % average of the 4 k's.
%         [k1 field(ti)] = feval(fderiv, ut, K, Vf, x, conjU(ti), Epenal)*dt;
%         k2 = feval(fderiv, ut + 0.5*k1, K, Vf, x, conjUhalf, Epenal)*dt;
%         k3 = feval(fderiv, ut + 0.5*k2, K, Vf, x, conjUhalf, Epenal)*dt;
%         k4 = feval(fderiv, ut + k3, K, Vf, x, conjU(ti + 1), Epenal)*dt;
%         U(:, ti + 1) = U(:, ti) + (k1 + 2*k2 + 2*k3 + k4)/6;
        t = t + dt;
    end
    tiIpln = (Nt-2):(Nt + 1);
    conjUhalf = NewtonIpln4(t0 + (tiIpln-1)*dt, conjU(:, tiIpln), t + 0.5*dt);
    [U(:, Nt + 1) field(Nt)] = unew(fderiv, U(:, Nt), K, Vf, x, conjU(:, Nt), conjUhalf, conjU(:, Nt + 1), Epenal, dt);
end

function [unext field] = unew(fderiv, ut, K, Vf, x, conjUlast, conjUhalf, conjUnext, Epenal, dt)
 % The k's are estimations of du. The result is composed of a weighted
 % average of the 4 k's.
    [deriv1 field] = fderiv(ut, K, Vf, x, conjUlast, Epenal);
    k1 = deriv1*dt;
    k2 = fderiv(ut + 0.5*k1, K, Vf, x, conjUhalf, Epenal)*dt;
    k3 = fderiv(ut + 0.5*k2, K, Vf, x, conjUhalf, Epenal)*dt;
    k4 = fderiv(ut + k3, K, Vf, x, conjUnext, Epenal)*dt;
    unext = ut + (k1 + 2*k2 + 2*k3 + k4)/6;
end