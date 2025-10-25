function u=matvec(v)
global A
global nvec
nvec=nvec+1;
u=A*v;
