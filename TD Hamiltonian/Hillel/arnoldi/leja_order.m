function z=leja_order(zset)
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
   if imag(z(j))~=0
       j=j+1;
       z(j,1)=conj(z(j-1,1));
       zv=abs(z(j)-zset).*zv;
   end
  [amax,imax]=max(zv);
  zv=zv/amax;
end