function u=matvec(v)
global A
global matvecs
matvecs=matvecs+1;
u=A*v;

