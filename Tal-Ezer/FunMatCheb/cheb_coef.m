function [c,R]=cheb_coef(fun,m,t,tol,a,b)
for i=1:m
    y(i)=cos((2*i-1)*pi/2/m);
    z(i)=.5*((b-a)*y(i)+a+b);
    f(i,1)=feval(fun,z(i)*t);
end
c=sqrt(2/m)*dct(f);
c(1)=c(1)/sqrt(2);
for j=1:m
    if abs(c(j)) < tol
        break
    end
end
c=c(1:j);
 