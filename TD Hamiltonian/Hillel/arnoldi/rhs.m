function fs=rhs(u,tvec,thalf)
global x
iu=sqrt(-1);
m=length(tvec);
for j=1:m
    fs(:,j)=-iu*((ff(tvec(j))-ff(thalf))*x).*u(:,j);
end
