% The program assumes the existence of the variable miu.
% Don't forget to change miu according to Vt3LS.
Vt = Vt3LS(100, 1000);
U = TDHexpH([1.1 2 3.1], Vt, [1;0;0], 100, 1000);
mmiu = zeros(1001, 1);
for ti=1:1001
    mmiu(ti) = exval(miu, U(:, ti));
end
mmiuw = dctI(mmiu);
Kre = 50;
dt = 0.5;
c = mmiu;
c = c(1:5:end);
FBDMrinp
%KBDMr