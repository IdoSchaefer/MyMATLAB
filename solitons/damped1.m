function DxDv = damped1(getxv, p)
x = getxv(1);
v = getxv(2);
g = p(1);
k = p(2);
Dx = v;
Dv = -g*v-k*x;
DxDv = [Dx; Dv];