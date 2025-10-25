function a= cheb(f)
[N,k]=size(f);   
for j=1:k
    a(:,j)=fft([f(N:-1:1,j); flipud(f(N-1:-1:2,j))]);               % Extend and compute fft
    a(1,j)=.5*a(1,j); a(N,j)=.5*a(N,j);
    a(:,j)=a(:,j)/(N-1);
end
a=a(1:N,:);
