function [Dchi field] = OCchim(chit, K, Vf, x, psit, Epenal)
    field = -imag(chit'*(x.*psit))/Epenal;
    Dchi = -1i*(ifft(K.*fft(chit)) + (Vf(x) - field*x).*chit);
end