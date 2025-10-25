function u=matvec(v)
global matvecs
global mu
matvecs=matvecs+1;
n2=length(v);
n=n2/2;
u(1:n,1)=v(n+1:n2,1);
u(n+1:n2,1)=fourdifft(v(1:n),2)-mu*v(n+1:n2,1);