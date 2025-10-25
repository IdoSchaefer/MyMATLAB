function [field, psi, relE] = OCTLS2(psi0, target, E, fguess, Epenal, T, dt, tol)
% fguess is a function handle.
    Nt = T/dt;
    field = zeros(1, Nt+1);
%    Vt = zeros(2, 2, Nt+1);
%    ezer = zeros(1, Nt+1);
    psi = zeros(2, Nt+1);
    chi = zeros(2, Nt+1);
    for ti = 1:(Nt+1)
        field(ti) = fguess((ti-1)*dt);
 %       Vt = [0,       -conj(field(ti));
 %            -field(ti), 0             ];
    end
    psi(:, 1) = psi0;
    lastfield = field;
    relE = tol + 1;
    maxNiter = 200;
    niter = 0;
    while relE>tol && niter<maxNiter
        sol = ode45(@TLSEf, [0 T], psi(:, 1), [], E, field, T);
        psi = deval(sol, 0:dt:T);
%        chi(:, Nt+1) = conj(target).*psi(:, Nt + 1);
        chi(:, Nt+1) = (target'*psi(:, Nt + 1))*target;
        sol = ode45(@TLSEf, [T 0], chi(:, Nt+1), [], E, field, T);
        chi(:, 1:Nt) = deval(sol, 0:dt:(T-dt));
        for ti = 1:(Nt+1)
            field(ti) = 1i/Epenal*(conj(chi(1, ti))*psi(2, ti) - conj(psi(1, ti))*chi(2, ti));
        end
%         integral = 0.5*conj(ezer(1))*ezer(1)*dt;
%         for ti = 2:Nt
%             integral = integral + conj(ezer(ti))*ezer(ti)*dt;
%         end
%         integral = integral + 0.5*conj(ezer(Nt+1))*ezer(Nt+1)*dt;
%         for ti = 1:(Nt+1)
%             field(ti) = ezer(ti)*sqrt(totalfE/integral);
% %            Vt(1, 2, ti) = -conj(field(ti));
% %            Vt(2, 1, ti) = -field(ti);
%         end
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
    end
    sol = ode45(@TLSEf, [0 T], psi(:, 1), [], E, field, T);
    psi = deval(sol, 0:dt:T);
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    niter
end