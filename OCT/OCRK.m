function [field, psi, relE, conv, niter] = OCRK(psi0, target, Vf, max_x, fguess, Epenal, T, dt, tol)
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
    field = zeros(1, Nt+1);
    lastfield = zeros(1, Nt+1);
    chi = zeros(Nx, Nt+1);
    psi = zeros(Nx, Nt+1);
tic
%     psisol = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
%    [t, psi] = ode45(@OCguess, 0:dt:T, psi0, [], K, Vf, x, fguess);
    [tpsi, psisol] = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
    psisol = psisol.';
toc
%    psiT = deval(psisol, T);
    psiT = psisol(:, end);
    for ti = 1:(Nt+1)
        lastfield(ti) = fguess((ti-1)*dt);
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
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        end
%        chisol = ode45(@OCchi, [T 0], chiT, [], K, Vf, x, psisol, Epenal);
%        [t, chi] = ode45(@OCchi, T:-dt:0, chiT, [], K, Vf, x, psi, Epenal, T, dt);
        [tchi, chisol] = ode45(@OCchi, [T 0], chiT, [], K, Vf, x, psisol, tpsi, Epenal);
        chisol = chisol.';
        chisol = chisol(:, end:-1:1);
        tchi = tchi(end:-1:1);
        for ti = 1:(Nt+1)
            chi(:, ti) = RK4interp(chisol, tchi, (ti - 1)*dt);
%             chi(:, ti) = RK4interp(chisol.y(:, end:-1:1), chisol.x(end:-1:1), (ti - 1)*dt);
        end
%        chi = deval(chisol, 0:dt:T);
        [tpsi, psisol] = ode45(@OCpsi, [0 T], psi0, [], K, Vf, x, chisol, tchi, Epenal);
        psisol = psisol.';
        for ti = 1:(Nt+1)
            psi(:, ti) = RK4interp(psisol, tpsi, (ti - 1)*dt);
        end
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
        for ti = 1:(Nt+1)
            field(ti) = -imag(chi(:, ti)'*(x.*psi(:, ti)))/Epenal;
        end
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
%    psi = deval(psisol, 0:dt:T);
    if nargin>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = target'*psiT;
        conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        conv = conv(:, 1:(niter + 1));
    end
    niter
end