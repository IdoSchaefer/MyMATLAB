function [T, R, total_norm, grad_norm, gradT, gradR] = TRcoef_temp(k0, V, xdomain, Nx)
% The function computes the transmission and reflection coefficients for a
% potential in a given spatial domain. It is assumed that the potential is
% 0 near both edges of the domain.
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
    if k0 == 0
        k0 = 10*eps;
    end
    E = k0^2/2;
    % Computation of the A(k_j) multiplication matrix B:
    B = eye(2);
    if nargout>3
        % The gradient matrices, wrt all V_j:
        gradB = zeros(2, 2, Nx);
    end
    for xi = 1:Nx
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
    % The relevant terms of the inverse transfer matrix:
    D11 = 0.5*exp(1i*k0*L)*(B(1, 1) + B(2, 2) + 1i*(B(1, 2)*k0 - B(2, 1)/k0));
    D21 = 0.5*exp(1i*k0*L)*(B(1, 1) - B(2, 2) + 1i*(B(1, 2)*k0 + B(2, 1)/k0));
    L_minus = D21/D11;
    R_plus = 1/D11;
    T = R_plus*conj(R_plus);
    R = L_minus*conj(L_minus);
    total_norm = T + R;
    if nargout>3
        % The relevant terms:
        gradD11(:, 1) = 0.5*exp(1i*k0*L)*(gradB(1, 1, :) + gradB(2, 2, :) + 1i*(gradB(1, 2, :)*k0 - gradB(2, 1, :)/k0));
        gradD21(:, 1) = 0.5*exp(1i*k0*L)*(gradB(1, 1, :) - gradB(2, 2, :) + 1i*(gradB(1, 2, :)*k0 + gradB(2, 1, :)/k0));
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