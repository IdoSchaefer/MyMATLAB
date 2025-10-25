function Dpsi = OCpsi(t, psit, K, Vf, x, chi, tchi, Epenal)
% function Dpsi = OCpsi(t, psit, K, Vf, x, chi, Epenal, T, dt)
%    chit = RK4interp(chisol.y(:, end:-1:1), chisol.x(end:-1:1), t);
%    T = 10; dt = 0.1;
    chit = RK4interp(chi, tchi, t);
    field = -imag(chit'*(x.*psit))/Epenal;
    Dpsi = -1i*(ifft(K.*fft(psit)) + (Vf(x) - field*x).*psit);
%            field];
end
