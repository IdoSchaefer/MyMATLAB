function u = chebpol(A, u0, T, m, leftb, rightb)
    A = A*T;
    leftb = T*leftb;
    rightb = T*rightb;
    y = cos(((1:m)*2 - 1)*pi/(2*m));
    z = 0.5*(y*(rightb - leftb) + rightb + leftb);
    c = dct(exp(z))/sqrt(m);
    c(2:m) = c(2:m)*sqrt(2);
%    c = chbcoef(@exp, leftb, rightb, m);
    v1 = u0;
    dim = size(A);
    I = eye(dim(1));
    Acheb = (2*A - (leftb + rightb)*I)/(rightb - leftb);
    v2 = Acheb*u0;
    u = c(1)*v1 + c(2)*v2;
    for k = 3:m
        vk = 2*Acheb*v2 - v1;
        u = u + c(k)*vk;
        v1 = v2;
        v2 = vk;
    end
end
        