x=-8:0.125:7.875;
X = diag(x);
k= -pi/0.125:2*pi/16:(pi/0.125 - 2*pi/16);
P=p2x(diag(k));
H = P^2/2 + X^2/2;
fi0 = 1/pi^(1/4)*exp(-x.^2/2)'/sqrt(8);
Vh = Vtx3(@(x, t) x*cos(t), x, 10, 100, 10);
U = TDHdiag_tscx(H, Vh, fi0, 10, 100, 10, 1e-2);
% Computation of the expectation value of x:
mx = zeros(1, 101);
for j=1:101
    mx(j) = exval(X, U(:, j));
end
t = 0:0.1:10;
figure
plot(t, mx)
error = mx - (-0.5*sin(t).*t);
figure
plot(t, error)
U1 = TDHdiag_tscx(H, Vh, fi0, 10, 100, 10, 1e-3);
U1(1:10, 1:5)