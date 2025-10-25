function [T, R, total_norm, grad_norm, gradT, gradR] = TRcoef2(k0, V, xdomain, Nx)
%%%%%%% Note: T and R are not the transmission and reflection coefficients,
%%%%%%% but the norms of the reflected and transmitted parts. It should be
%%%%%%% fixed.
%%%% Note: the gradients are computed correctly only for real final k. It
%%%% should be fixed.
%%% I looked at it after more than a year, and I think that this note is
%%% incorrect. Everything looks fine.
% The function computes the transmission and reflection coefficients for a
% potential in a given spatial domain. It is assumed that the potential is
% 0 near the left edge of the domain.
% Input:
% k0: The wave number of the incident wave.
% V: A function handle; the potential function. It is possible to insert
% the potential value vector itself instead.
% xdomain: The domain of the x grid
% Nx: The number of eqidistant points in the x grid. Needed to define the
% accuracy of the computation.
% Output:
% T: The transmission coefficient
% R: The reflection coefficient
% total_norm: The total norm of the reflected and transmitted waves. Should
% be 1 for a real potential.
% grad_norm: The gradient of the norm with respect to the V values. Should
% be 0 for a real potential. The 1'st Nx values are the gradient wrt the
% real part of the potential, and the last Nx values are wrt the imaginary
% part.
% gradT: The gradient of T with respect to the V values.
% gradR: The gradient of R with respect to the V values.
    L = xdomain(2) - xdomain(1);
    dx = L/(Nx - 1);
    x = (xdomain(1):dx:xdomain(2)).';
    if length(V) == 1
        % If V is a function handle:
        V = V(x);
    end
    E = k0^2/2;
    % Computation of the A(k_j) multiplication matrix B:
    B = eye(2);
    if nargout>3
        % The gradient matrices, wrt all V_j (excluding V_N):
        gradB = zeros(2, 2, Nx - 1);
    end
    for xi = 1:(Nx - 1)
        k = sqrt(2*(E - V(xi)));
        if k~=0
            A = [cos(k*dx)    -sin(k*dx)/k; %-dx*sincm(k*dx); Not recommended. It takes more than half of the computational time.
                 k*sin(k*dx)   cos(k*dx)  ];            
        else
            A = [1  -dx; 
                 0   1  ];
        end
        if nargout>3
            for ki = 1:(xi - 1)
                gradB(:, :, ki) = gradB(:, :, ki)*A;
            end
% In the earlier version, the gradient was computed wrt k. I don't know
% exactly why, but this gives a wrong picture on the imaginary k values. 
%             gradB(:, :, xi) = B*[-dx*sin(k*dx),                 (sin(k*dx) - k*dx*cos(k*dx))/k^2;
%                                  sin(k*dx) + k*dx*cos(k*dx),    -dx*sin(k*dx)];
% In the current version, the gradient is wrt V. This ammounts of
% multiplication by dk/dV = -1/k.
            if k~=0
                gradB(:, :, xi) = B*[dx*sin(k*dx)/k,                 (k*dx*cos(k*dx) - sin(k*dx))/k^3;
                                     -sin(k*dx)/k - dx*cos(k*dx),    dx*sin(k*dx)/k                     ];
            else
                gradB(:, :, xi) = B*[dx^2,     -dx^3/3;
                                     -2*dx,    dx^2    ];
            end                                                
        end
        B = B*A;
    end
    if E~=V(Nx)
        kN = sqrt(2*(E - V(Nx)));
    else
        kN = sqrt(2*(10*eps));
        %kN = 1e5*eps;
    end
    % The relevant terms of the inverse transfer matrix:
    exp_ikNL_d2k0 = exp(1i*kN*L)/(2*k0);
    D11 = exp_ikNL_d2k0*(k0*B(1, 1) + kN*B(2, 2) + 1i*(k0*kN*B(1, 2) - B(2, 1)));
    D21 = exp_ikNL_d2k0*(k0*B(1, 1) - kN*B(2, 2) + 1i*(k0*kN*B(1, 2) + B(2, 1)));
    L_minus = D21/D11;
    R_plus = exp(1i*kN*L)/D11;
    T = R_plus*conj(R_plus);
    R = L_minus*conj(L_minus);
    total_norm = T + R;
    if nargout>3
        % The relevant terms:
        gradD11 = zeros(Nx, 1);
        gradD21 = zeros(Nx, 1);
        gradD11(1:(Nx - 1)) = exp_ikNL_d2k0*(k0*gradB(1, 1, :) + kN*gradB(2, 2, :) + 1i*(k0*kN*gradB(1, 2, :) - gradB(2, 1, :)));
        gradD11(Nx) = -exp_ikNL_d2k0*(B(2, 2) - L*(k0*kN*B(1, 2) - B(2, 1)) + 1i*(L*(k0*B(1, 1) + kN*B(2, 2)) + k0*B(1, 2)))/kN;
        gradD21(1:(Nx - 1)) = exp_ikNL_d2k0*(k0*gradB(1, 1, :) - kN*gradB(2, 2, :) + 1i*(k0*kN*gradB(1, 2, :) + gradB(2, 1, :)));
        gradD21(Nx) = -exp_ikNL_d2k0*(-B(2, 2) - L*(k0*kN*B(1, 2) + B(2, 1)) + 1i*(L*(k0*B(1, 1) - kN*B(2, 2)) + k0*B(1, 2)))/kN;
        gradL_minus = (gradD21*D11 - D21*gradD11)/D11^2;
        %gradL_minus = gradD21/D11 - D21*gradD11/D11^2;
        gradR_plus = -gradD11/D11^2;
        conjRp_gradRp = conj(R_plus)*gradR_plus;
        conjLm_gradLm = conj(L_minus)*gradL_minus;
        gradT = 2*[real(conjRp_gradRp); -imag(conjRp_gradRp)];
        gradR = 2*[real(conjLm_gradLm); -imag(conjLm_gradLm)];
%         gradT = 2*real(conj(R_plus)*gradR_plus);
%         gradR = 2*real(conj(L_minus)*gradL_minus);
        grad_norm = gradT + gradR;
    end
end