%This program assumes the existance of a variable: mind, which is the
%value of the minimal amplitude of the terms to be saved and displayed.
mind = input('Enter the minimal amplitude: ');
d_fil = d_result(d_result > mind);
w_fil = w_result(d_result > mind);
%[w_fil ow_fil] = sort(w_fil);
%wn_fil = real(w_fil)*1e15/(3e10*2*pi)
%[wn_fil, ow_fil] = sort(real(w_fil)*1e15/(3e10*2*pi))
wn_fil = real(w_fil)*1e15/(3e10*2*pi)
gamma_fil = imag(w_fil)
d_fil