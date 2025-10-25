function v=chbeig(R,t,matvec,n,m,s)
for j=1:m
    y(j)=cos((2*j-1)*pi/(2*m));
    z(j)=R*(y(j)-1)/2;
end
for i=1:s
    tt(i)=i*t/s;
    for j=1:m
        f(j,i)=exp(R*tt(i)*z(j));
    end
    c(:,i)=sqrt(2/m)*dct(f(:,i));
    c(1,i)=c(1,i)/sqrt(2);
end
v1=rand(n,1); 
v2=(2/R)*feval(matvec,v1)+v1;
for i=1:s
    v(:,i)=c(1,i)*v1+c(2,i)*v2;
end
for k=3:m
    v3=(2/R)*feval(matvec,v2)+v2;
    v3=2*v3-v1;
    v1=v2;
    v2=v3;
    for i=1:s
        v(:,i)=v(:,i)+c(k,i)*v3;
    end
end
