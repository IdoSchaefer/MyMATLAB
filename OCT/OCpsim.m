function [Dpsi field] = OCpsim(psit, K, Vf, x, chit, Epenal)
    field = -imag(chit'*(x.*psit))/Epenal;
    Dpsi = -1i*(ifft(K.*fft(psit)) + (Vf(x) - field*x).*psit);
end