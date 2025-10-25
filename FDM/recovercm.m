%This program assumes the existance of a variable: mind, which is the
%value of the minimal amplitude of the terms to be saved and displayed.
mind = input('Enter the minimal amplitude: ');
meanexpd = abs(expd).*1./(imag(expw)*tmax).*(exp(imag(expw)*tmax) - 1);
expd_fil = expd(meanexpd > mind);
meand_fil = meanexpd(meanexpd > mind);
expw_fil = expw(meanexpd > mind);
% optional - flipping of positive "damping" terms:
%expw_fil(imag(expw_fil)>0) = conj(expw_fil(imag(expw_fil)>0));
% another option:
%expw_fil(imag(expw_fil)>0) = real(expw_fil(imag(expw_fil)>0));
[meand_fil, omd_fil] = sort(meand_fil);
expw_fil = expw_fil(omd_fil);
expd_fil = expd_fil(omd_fil);
wn_fil = real(expw_fil)*1e15/(3e10*2*pi)
% optional: to display the damping coefficients:
gamma_fil = imag(expw_fil)
meand_fil
%crec = zeros(2*K, 1);
%for n = 0:(2*K - 1)
    % The time index is n+1.
%    crec(n+1) = sum(expd_fil.*exp(-i*expw_fil*dt*n));
%end
% The same thing, a better option:
%n = 0:2*K-1;
%crec = expd_fil*exp(-i*expw_fil.'*dt*n);
n = 0:(2*K-1)*dt;
crec = expd_fil*exp(-i*expw_fil.'*n);