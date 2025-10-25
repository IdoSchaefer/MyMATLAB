function result = Hpsicheb(K, V, psi, leftb, rightb)
% The operation of the Hamiltonian in the domain [-1 1], for the use of Chebychev expansion.
% Look in Hpsi for more details.
% The original domain of the eigenvalues is: [leftb rightb].
    result = (2*(ifft(K.*fft(psi)) + V.*psi) - (leftb + rightb)*psi)/(rightb - leftb);
end