syms x y z a b c d;
f=a*x*0.5*(1-tanh(b*(x-c)^d));
Df=diff(f);
D2f=diff(Df);
D3f=diff(D2f);
D4f=diff(D3f);
Df0=subs(Df,x,0);
D2f0=subs(D2f,x,0);
D3f0=subs(D3f,x,0);
D4f0=subs(D4f,x,0);
params=solve(Df0-0.1926, D2f0-0.01763, D3f0+0.2233, D4f0+0.2870);
load HCl
miusymf=subs(f, {a,b,c,d}, {params.a(3) params.b(3) params.c(3) params.d(3)});
miuv=subs(miusymf, x);
miuv=real(miuv);
miuf1 = @(x) 0.1926*x + 1/2*0.01763*x.^2 - 1/6*0.2233*x.^3 -1/24*0.287*x.^4;
ap = params.a(3);
bp=params.b(3);
cp=params.c(3);
dp=params.d(3);
ap=real(ap); bp=real(bp);
ap=subs(ap);bp=subs(bp);cp=subs(cp);dp=subs(dp);
miuf2=@(x) ap*x*0.5.*(1-tanh(bp*(x-cp).^dp));
% miuf2 is different from miuv. The 2 functions differ only in the region
% in which miuf1 is also very different.
bpc = params.b(3);
bpc=subs(bpc);
miuf3=@(x) ap*x*0.5.*real(1-tanh(bpc*(x-cp).^dp));
% miuf3 represents miuv.