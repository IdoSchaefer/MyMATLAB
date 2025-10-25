[B, Du] = eig(U1, U0);
w_result = -log(diag(Du).')/(i*dt);
w_result(real(w_result)<0) = w_result(real(w_result)<0) + 2*pi/dt;
[w_result, orderw] = sort(w_result);
rwerr = real(w_result - w)./real(w);
imwerr = imag(w_result - w)./imag(w);
werr = sqrt(rwerr.^2 + imwerr.^2);
