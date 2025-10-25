function Phi_out = multigrid(omega, h, Phi_in, fv, Niter, min_level)
    N = sqrt(length(Phi_in)) + 1;
    level = log2(N);
    L = Laplacian_matrix(N, h);
    if level == min_level
        Phi_out = L\fv;
        %Phi_out = smooth(update_matrix(N, omega), L, omega, h, Phi_in, fv, 10*Niter);
        return
    end
    Phi_out = smooth(update_matrix(N, omega), L, omega, h, Phi_in, fv, Niter);
    r = fine2coarse(fv - L*Phi_out);
    Ncoarse = N/2;
    psi = multigrid(omega, h, zeros((Ncoarse - 1)^2, 1), r, Niter, min_level);
    Phi_out = Phi_out + coarse2fine(psi);
    Phi_out = smooth(update_matrix(N, omega), L, omega, h, Phi_out, fv, Niter);
end