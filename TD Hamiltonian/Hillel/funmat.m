function u=funmat(v,A,t,fun)
[vec,eg]=eig(A);
eg=diag(eg);
n=length(v);
for j=1:n
    w(j,1)=feval(fun,eg(j),t);
end
u=vec*(diag(w))*(vec\v);