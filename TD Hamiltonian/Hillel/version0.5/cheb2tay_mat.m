function A=cheb2tay_mat(m,dt)
u=zeros(m,1);
for j=1:m
    u(j)=1;
    A(j,:)=cheb2tay(u,dt);
    u(j)=0;
end
    