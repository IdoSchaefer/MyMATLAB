function u=matvec(v)
global x
global thalf
global matvecs
matvecs=matvecs+1;
iu=sqrt(-1);
fs=ff(thalf)*(x.*v);
u=(iu/2)*(fourdifft(v,2)-(x.*x).*v);
u=u-iu*fs;

