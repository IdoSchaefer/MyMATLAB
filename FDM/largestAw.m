% (The program assumes the existence of a variable:) Ndisp, which is the
% number of the frequencies with the largest amplitudes to be displayed.
Ndisp = input('Enter the number of the frequencies with maximal amplitude: ');
[d_resulto, orderd] = sort(d_result);
%wn_od = wave_number(orderd);
w_od = wreal(orderd);
Nd = length(d_result);
%wndisp = wn_od((Nd-Ndisp+1):Nd)
wdisp = w_od((Nd-Ndisp+1):Nd)
ddisp = d_resulto((Nd-Ndisp+1):Nd)