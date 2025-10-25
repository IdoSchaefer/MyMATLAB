function a= cheb_row(f)
[N,k]=size(f);   
for j=1:k
    a(j,:)=fft([f(N:-1:1,j); flipud(f(N-1:-1:2,j))]);               % Extend and compute fft
    a(j,1)=.5*a(j,1); a(j,N)=.5*a(j,N);
    a(j,:)=a(j,:)/(N-1);
end
a=a(:,1:N);
