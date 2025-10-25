function [r,normr,nre]=reorth(q,r,normr,index,alpha,method)
[n,k1]=size(q);
k=length(index);
if k==k1 & index(:)==[1:k]';
    simple=1;
else
    simple=0;
end
normr_old=normr;
if method==1
    if simple 
        r=r-q*(q'*r);
    else
        r=r-q(:,index)*(q(:,index)'*r);
    end
else
    for i=index
        r=r-q(:,i)*(q(:,i)'*r);
    end
end
nre=1;
normr=norm(r);
if normr < alpha*normr_old
    if method==1
        if simple
            r=r-q*(q'*r);
        else
            r=r-q(:,index)*(q(:,index)'*r);
        end
    else
        for i=index
            r=r-q(:,i)*(q(:,i)'*r);
        end
    end
    normr_old=normr
    normr=norm(r);
    if normr < alpha*normr_old
        %r is in span(q) to full accuracy hence accept r = 0 as the new
        %vector
        r=zeros(n,1);
        normr=0;
    end
    nre=2;
end