function u=funmat(fun,A,t,w)
[vec,d]=eig(full(A));
d=diag(d);
n=length(d);
for i=1:n
    fd(i,1)=feval(fun,d(i),t);
end
ww=vec\w;
u=vec*(fd.*ww);