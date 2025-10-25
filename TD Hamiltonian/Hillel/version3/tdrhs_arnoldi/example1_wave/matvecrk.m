function u=matvecrk(v)
global matvecs
global x
matvecs=matvecs+1;
iu=sqrt(-1);
u=(iu/2)*(fourdifft(v,2)-(x.*x).*v);
