function c=chebcoef(a,b,m1,m2,fun,tvec)
for i=1:m2
     y=-cos((i-1)*pi/(m2-1));
     z(i)=(1/2)*((b-a)*y+a+b);
     for k=1:m1
         f(i,k)=feval(fun,z(i),tvec(k),m1);
     end
end
c=cheb(f);
