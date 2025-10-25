function DxDv = harmonicf(getxv, sqromega)
x = getxv(1);
v = getxv(2);
Dx = v;
Dv = -sqromega*x;
DxDv = [Dx; Dv];