%This program assumes the existance of a variable: minw, which is the
%value of the minimal angular freqency to be saved and displayed.
w_fil = w_result(w_result > minw);
d_fil = d_result(w_result > minw);
[w_fil ow_fil] = sort(w_fil);
w_fil
d_fil = d_fil(ow_fil);