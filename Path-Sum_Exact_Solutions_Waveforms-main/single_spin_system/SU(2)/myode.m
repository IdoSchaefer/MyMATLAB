function dg = myode(t,para,g)

[Cx,Cy] = chirp_fun(t,para);

H = [para.Omega/2,   Cx-1i*Cy;...
     Cx+1i*Cy,      -para.Omega/2];

H_hat = kron(eye(2),H)-kron(conj(H),eye(2));

%gvec = [g(1); g(2); g(3); g(4)];

dg = -1i*H_hat*g;

end