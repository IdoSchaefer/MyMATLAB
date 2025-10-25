function u=matvec(v)
global A
global iter
iter=iter+1;
u=A*v;
