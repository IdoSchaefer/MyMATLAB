function result = chebc2result(Ccheb, xdomain, xresult)
    xrcheb = (2*xresult - xdomain(1) - xdomain(2))/(xdomain(2) - xdomain(1));
    m = length(Ccheb);
    t1 = ones(size(xresult));
    t2 = xrcheb;
    result = Ccheb(1)*t1 + Ccheb(2)*t2;
    for k = 3:m
          tk = 2*xrcheb.*t2 - t1;
          result = result + Ccheb(k).*tk;
          t1 = t2;
          t2 = tk;
    end
end