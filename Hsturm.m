function H = Hsturm(x, p, Vf, b, m)
% The function computes the "Hamiltonian" for a Sturmian basis.
    if nargin < 5
        m = 1;
    end
    Nx = length(x);
    Kp = p.^2/(2*m);
    bvec = ones(Nx, 1)*b;
%    bvec(1:floor(Nx/2)) = -bvec(1:floor(Nx/2));
    bvec(x<0) = -bvec(x<0);
    H = diag(Vf(x) - b^2/(2*m)) + Nx*(ifft(ifft(diag(Kp))')' + 1i*diag(bvec)*ifft(ifft(diag(p/m))')');
end