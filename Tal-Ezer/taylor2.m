function u = taylor2(A, u0, T, Ndt)
    dt = T/Ndt;
    dim = size(A);
    I = eye(dim(1));
    u = (I + dt*A + 0.5*(dt*A)^2)^Ndt*u0;
end