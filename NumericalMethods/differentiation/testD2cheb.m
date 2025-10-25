function [cD2, D2fun, c] = testD2cheb(f, leftb, rightb, N)
    c = chebc(f, leftb, rightb, N);
    cD2 = D2chebc(c, rightb - leftb);
    D2fun = idct([cD2(1), cD2(2:N)/sqrt(2)]*sqrt(N));
end