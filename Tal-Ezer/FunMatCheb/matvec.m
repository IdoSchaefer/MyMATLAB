function u=matvec(v)
global x g gtag
global nvec
nvec=nvec+1;
n=length(v);
% w=drv(v);
% u=drv(w);
 u=drv2(v);
% u=sqrt(-1)*(u+(x.^2).*v);
u=sqrt(-1)*(u/2-(x.^2)/2.*v);
