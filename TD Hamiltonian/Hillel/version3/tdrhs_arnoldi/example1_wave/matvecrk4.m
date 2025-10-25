function u=matvecrk4(v)
global matvecs
global dx
matvecs=matvecs+1;
n=length(v);
for j=2:n-1
    u(j,1)=(v(j+1)-v(j-1))/2/dx;
end
u(1)=0;
u(n)=(v(n)-v(n-1))/dx;