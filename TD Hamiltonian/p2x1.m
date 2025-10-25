function Mx = p2x1(Mp)
% The function gets a matrix in the p space (in atomic units), and returns the corresponding  matrix in the x space.
% The p domain is assumed to be: [-Pmax, Pmax), and the matrix Mp is ordered 
% in the appropriate way for a Fourier transform without shifting.
    szMp = size(Mp);
    N = szMp(1);
    Mx = N*ifft(ifft(Mp)')';
end