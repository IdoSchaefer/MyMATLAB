function fs=rhs0(tvec,thalf)
iu=sqrt(-1);
m=length(tvec);
for j=1:m
    fs(j)=-iu*(ff(tvec(j))-ff(thalf));
end
