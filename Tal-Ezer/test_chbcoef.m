function test_chbcoef
m=40;
a=-100; b=0;
c=chbcoef('exp',a,b,m)

function c=chbcoef(fun,a,b,m)
for j=1:m
    y=cos((2*j-1)*pi/(2*m));
    z(j)=.5*((b-a)*y+a+b);
    f(j,1)=feval(fun,z(j));
  %  f(j,1)=feval(fun,y);
end
c=sqrt(2/m)*dct(f);
c(1)=c(1)/sqrt(2);

