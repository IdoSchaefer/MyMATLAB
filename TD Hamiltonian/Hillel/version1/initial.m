function [x,u0,L]=initial(n)
xmin=-10;
dx=20/n;
L=abs(xmin);
for i=1:n
    x(i,1)=xmin+(i-1)*dx;
    u0(i,1)=exp(-x(i)^2);
end
