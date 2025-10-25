function A=matrix(n)
w=ones(n,1);z=zeros(n,1);
z=-[1:n]';z=z/n;z=zeros(n,1);
A=spdiags([ w -2*w w ],[-1:1],n,n);
A=spdiags([ -w z w ],[-1:1],n,n);

