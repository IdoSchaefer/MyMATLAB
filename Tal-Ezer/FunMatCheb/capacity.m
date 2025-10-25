function  ro=capacity(znew,z)
n=length(z);
m=length(znew);
ro=1;
if n ==0
    ro=max(abs(znew))/4;
    return
end
r=ones(m,1);
for k=1:n
    r=abs(znew-z(k)).^(1/n).*r;
end
ro=max(r);
