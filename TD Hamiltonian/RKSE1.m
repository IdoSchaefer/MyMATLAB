function U = RKSE1(H0, Vt, u0, T, dt)
%%% not written yet.
% The program solves The Schrodinger equation, with a Hamiltonian of the
% form: H = H0 + Vt.
% H0 is a matrix.
% Vt is a 3d tensor, when the 3'rd index represent the time index, and the first 2 indices are the matrix. 
% u0 is the initial state vector.
% T is the total time.
% dt is the time difference between the time points that will be returned.
% The columns of the solution, U, represent different time points.
    if nargin<5
        dt = 0.1;
    end
    szx = size(x);
    if szx(2)>1
        x = x.';
    end
    sol = ode45(@SE, [0 T], u0, [], H0, Vt, x);
    U = deval(sol, 0:dt:T);
end