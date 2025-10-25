function u=new_guess(v,w,m1,m2,dt,tvec,a,b,fun)
tj=ones(m1,1);
for j=1:m1-1
    tj(:,j+1)=tj(:,j).*(dt+tvec);
end
c=chebcoef(a,b,m1,m2,fun,dt+tvec);
u=v(:,1:m1)*tj'+w*c;
