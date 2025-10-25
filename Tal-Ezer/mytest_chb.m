function [result, relE] = mytest_chb(fun, x, a, b, m)
    c = chbcoef(fun, a, b, m);
    y = (2*x - a - b)/(b-a);
    t1 = 1; t2=y;
    result= c(1) + c(2)*t2;
    for k=3:m
          tk=2*y*t2 - t1;
          result =result+c(k)*tk;
          t1=t2; t2=tk;
    end
    ex_result = feval(fun, x);
    relE = (result - ex_result)/ex_result;

function c=chbcoef(fun,a,b,m)
for j=1:m
    y=cos((2*j-1)*pi/(2*m));
    z(j)=.5*((b-a)*y+a+b);
    f(j,1)=feval(fun,z(j));
end
c=sqrt(2/m)*dct(f);
c(1)=c(1)/sqrt(2);