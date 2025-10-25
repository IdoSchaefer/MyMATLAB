function b=convrt(z,a)
[n,m]=size(a);
b=[];
c(1)=-z(1); c(2)=1;
b(:,1)=a(:,1)-a(:,2)*z(1); b(:,2)=a(:,2);
for j=2:m-1
    d(1)=-z(j)*c(1);
    for k=2:j
        d(k)=-z(j)*c(k)+c(k-1);
    end
    d(j+1)=c(j);
    for k=1:j
        b(:,k)=b(:,k)+a(:,j+1)*d(k);
    end
    b(:,j+1)=a(:,j+1)*d(j+1);
    c=d;
end

