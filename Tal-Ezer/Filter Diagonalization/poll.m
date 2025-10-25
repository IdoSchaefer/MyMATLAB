function z=poll(y,c,R)
m=length(c);
t1=1; 
t2=2*y/R+1;
z=c(1)*t1+c(2)*t2;
for k=3:m
    t3=2*(2*y/R+1)*t2-t1;
    t1=t2;
    t2=t3;
    z=z+c(k)*t3;
end
    