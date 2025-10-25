function fs=rhs(u,tvec,thalf)
global x
iu=sqrt(-1);
[n,m]=size(u);
xs=x.*x;
for j=1:m
    fs(:,j)=-iu*((ff(tvec(j))-ff(thalf))*x).*u(:,j);
end
