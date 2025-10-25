function relE = Echeb(f, x, leftb, rightb, N)
    c = chebc(f, leftb, rightb, N);
    xcheb = (2*x - leftb - rightb)/(rightb-leftb);
    t1 = 1; t2=xcheb;
    result= c(1) + c(2)*t2;
    for k=3:N
          tk=2*xcheb*t2 - t1;
          result =result+c(k)*tk;
          t1=t2; t2=tk;
    end
    ex_result = feval(f, x);
    relE = (result - ex_result)/ex_result;
end