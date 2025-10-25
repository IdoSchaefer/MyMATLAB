function [field, psi, conv, niter, mallniterc] = OCBFGSMLS(psi0, target, E0, Edomain, miu, fguess, Epenal, T, dt, Nt_ts, Ncheb, tol)    
    Nt = T/dt;
    dim = length(psi0);
    maxNiter = 500;
    tolprop = 1e-2*tol;
    conv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
    summniterc = 0;
    nprop = 0;
    field = zeros(1, Nt + 1);
    psi = zeros(dim, Nt + 1);
    allfield = zeros(1, Nt*(Nt_ts - 1) + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    for tsi = 1:Nt
        for ti = 1:(Nt_ts - 1)
            allfield((tsi - 1)*(Nt_ts - 1) + ti) = fguess((tsi - 1)*dt + t_ts(ti));
        end
    end
    allfield(Nt*(Nt_ts - 1) + 1) = fguess(Nt*dt);
    integw = chebgridw(Nt, Nt_ts, dt);
%    options = optimset('gradObj', 'on', 'OutputFcn', @getconv, 'TolX', tol, 'MaxIter', maxNiter);  
%    options = optimset('gradObj', 'on', 'OutputFcn', @getconv, 'MaxIter', maxNiter);
    % This line was changed in 24/12/17, more than six years after the
    % original one:
    options = optimoptions('fminunc', 'algorithm', 'quasi-newton', 'OptimalityTolerance', tol, 'StepTolerance', tol, 'SpecifyObjectiveGradient',true,...
        'OutputFcn', @getconv, 'MaxIter', maxNiter);
    [allfield, lastJval, exitflag, output] = fminunc(@Jeval, allfield, options);
    %exitflag
    output
    niter = output.iterations;
    conv = conv(1:niter+1);
    [allpsi, field(1:Nt)] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    field(Nt + 1) = allfield(Nt*(Nt_ts - 1) + 1);
    mallniterc = summniterc/nprop;
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad] = Jeval(allfield)
        [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        overlap = target'*allpsi(:, Nt_ts, Nt);
        chiT = overlap*target;
        [allchi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield((Nt*(Nt_ts - 1) + 1):-1:1));
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        minusJ = -overlap*conj(overlap) + Epenal*sum(integw.*allfield.^2);
% The integration can be performed more efficiently, using the program:
% chebweigts. I've no time at the moment.
%         for tsi = 1:Nt
%             minusJ = minusJ + Epenal*intchebpb(allfield((tsi - 1)*(Nt_ts - 1) + (1:Nt_ts)).^2, dt);
%         end
        minusgrad = zeros(1, Nt*(Nt_ts - 1) + 1);
        for tsi = 1:Nt
            for ti = 1:(Nt_ts - 1)
                minusgrad((tsi - 1)*(Nt_ts - 1) + ti) = 2*(allfield((tsi - 1)*(Nt_ts - 1) + ti)*Epenal +...
                    imag(allchi(:, Nt_ts - ti + 1, Nt - tsi + 1)'*miu*allpsi(:, ti, tsi)));
            end
        end
        minusgrad(Nt*(Nt_ts - 1) + 1) = 2*(allfield(Nt*(Nt_ts - 1) + 1)*Epenal +...
                imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt)));
    end
    
    function stop = getconv(allfield, optimValues, state)
        stop = false;
        optimValues.iteration
        conv(optimValues.iteration + 1) = -optimValues.fval;
    end
end