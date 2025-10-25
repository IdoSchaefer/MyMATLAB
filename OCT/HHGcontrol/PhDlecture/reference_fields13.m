[fieldtcw1b, fieldwcw1b, psicw1b, evatcw1b, evawcw1b, evmiutcw1b, evmiuwcw1b, mniterccw1b, Jcw1b, J1cw1b, J2cw1b, ~, Jpnormcw1b] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*max(abs(fieldt1b))/max(abs(cw_con)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[fieldtcwe1b, fieldwcwe1b, psicwe1b, evatcwe1b, evawcwe1b, evmiutcwe1b, evmiuwcwe1b, mniterccwe1b, Jcwe1b, J1cwe1b, J2cwe1b, ~, Jpnormcwe1b] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*sqrt(sum(fieldt1b.^2)/sum(cw_con.^2)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
%%% Ionization:
%%% sqnorm(psicw1b(:,end)): 7.8340e-01
%%% sqnorm(psicwe1b(:,end)): 9.9879e-01
%%% Fluence:
%%% sum(fieldt1b.^2)*0.2: 9.9217e-01
%%% sum(fieldtcw1b.^2)*0.2: 2.6544e+00

