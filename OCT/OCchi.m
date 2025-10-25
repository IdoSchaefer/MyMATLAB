function Dchi = OCchi(t, chit, K, Vf, x, psi, tpsi, Epenal)
% function Dchi = OCchi(t, chit, K, Vf, x, psi, Epenal, T, dt)
    psit = RK4interp(psi, tpsi, t);
%    psit = RK4interp(psisol.y, psisol.x, t);
    field = -imag(chit'*(x.*psit))/Epenal;
    Dchi = -1i*(ifft(K.*fft(chit)) + (Vf(x) - field*x).*chit);
end
