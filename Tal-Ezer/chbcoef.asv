function c=chbcoef(fun,a,b,m)
for j=1:m
    y=cos((2*j-1)*pi/(2*m));
    z(j)=.5*((b-a)*y+a+b);
    f(j,1)=feval(fun,z(j));
end
c=sqrt(2/m)*dct(f);
c(1)=c(1)/sqrt(2);