function A=problem1(n)
w=ones(n,1);
A=spdiags([-w 2*w -w],-1:1,n,n);

