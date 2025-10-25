function test_gs
n=1000;m=40;
u(:,1)=rand(n,1);
eps=1.e-8;
for j=2:m
    u(:,j)=u(:,1)+eps*rand(n,1);
end
[v,r]=gs(u);
norm(eye(m)-v'*v)
[v,r]=mgs(u);
norm(eye(m)-v'*v)
