function z=leja_order(zset)
[k1,k2]=size(zset);
if k1==1
    zset=zset';
end
m=length(zset);
mmj=m;
zv=ones(mmj,1);
[a1,imax]=max(abs(zset));
j=0;
while j < m && mmj > 0
   j=j+1;
   z(j,1)=zset(imax);
   zset=zset([1:imax-1 imax+1:mmj]);
   zv=zv([1:imax-1 imax+1:mmj]);
   mmj=mmj-1;
   if mmj==0
       return
   end
   zv=abs(z(j)-zset).*zv;
   [amax,imax]=max(zv);
end
