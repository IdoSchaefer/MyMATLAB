function fs=rhs1(u,tvec)
global t
global H
% global mu
% global E0
iu=sqrt(-1);
m=length(tvec);
for j=1:m
    HT = diag([1 2])*(-1i);
    HT(1, 2) = -1i*exp(1i*(t+tvec(j)));
    HT(2, 1) = -1i*exp(-1i*(t+tvec(j)));    
%     HT=zeros(2);
%     HT(1,2)=-iu*mu/2*E0*S(t+tvec(j));
%     HT(2,1)=-iu*mu/2*E0*S(t+tvec(j));
    fs(:,j)=(HT-H)*u(:,j);
end
 