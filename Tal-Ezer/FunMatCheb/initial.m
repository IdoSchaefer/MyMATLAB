function u=initial(n)
ip=sqrt(-1);
for i=1:n
    %x(i,1)=-pi+(i-1)*2*pi/n;
    x(i,1)=-8+(i-1)*2*8/n;
    u(i,1)=exp(ip*x(i))*pi^(-1/4)*exp(-x(i)^2/2);
    abs(u(i));
end
ip