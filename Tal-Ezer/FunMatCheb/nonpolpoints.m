function [x,g,gtag]=nonpolpoints(p,n)
% The comments are mine:
for j=1:n
    y(j,1)=cos((j-1)*pi/(n-1));
    x(j,1)=(1/p)*asin(y(j)*sin(p));
    % g is dy/dx. The two ways are equivalent. One is 1/(dx/dy) in the
    % terms of y:
%      g(j,1)=p/sin(p)*sqrt(1-(y(j)*sin(p))^2);
% The other way is dy/dx in the terms of x:
     g(j,1)=p/sin(p)*cos(p*x(j));
     % This is (d^2)y/(dx)^2:
     gtag(j,1)=-p^2/sin(p)*sin(p*x(j));
end
