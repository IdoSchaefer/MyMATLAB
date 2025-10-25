function fs=rhs(u,un,tvec)
global x
global mu
m=length(tvec);
n=length(u);
fs=zeros(n,m);
for j=1:m
    fs(n/2+1:n,j)=mu*10*cos(10*(x+tvec(j)));
end
