T = 1e3;
[allNt5, allmv5, aller5, max_ers5] = errorSGarticle(T, 5, 5, 4e3, 17);
pause(0.1)
[allNt7, allmv7, aller7, max_ers7] = errorSGarticle(T, 7, 7, 2.85e3, 14);
pause(0.1)
[allNt9, allmv9, aller9, max_ers9] = errorSGarticle(T, 9, 9, 2.85e3, 11);
% Slopes:
% polyfit(log10(allmv5(3:12)), log10(aller5(3:12)), 1)
% ans = -8.9783e+00   3.9600e+01
% 
% polyfit(log10(allmv7(3:12)), log10(aller7(3:12)), 1)
% ans = -8.7738e+00   3.7185e+01
% 
% polyfit(log10(allmv9(3:10)), log10(aller9(3:10)), 1)
% ans =  -1.0579e+01   4.6369e+01
