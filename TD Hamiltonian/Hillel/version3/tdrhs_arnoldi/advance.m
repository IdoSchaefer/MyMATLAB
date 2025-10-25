function   u=advance(v,matvec,m1,m2,tvec,s,tol,cycles,ireal)
tj=ones(m1,1);
u=v*tj';
for j=1:m1-1
    tj=tj.*tvec;
    v=(feval(matvec,v)+s(:,j))/j;
    u=u+v*tj';
end
v=(feval(matvec,v)+s(:,m1))/m1;
w=newton(v,tvec',matvec,'fun',m2,tol,cycles,ireal);
u=u+w;