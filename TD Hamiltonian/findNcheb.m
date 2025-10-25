function Ncheb = findNcheb(f, leftb, rightb, OE, Nstart)
% The function finds the number of Chebyshev coefficients that have a
% greater value than OE, the order of the desired error.
% Nstart is an optional input argument, for the number of coefficients to start from;
% Use it if the program fails to give a resonable result.
    if nargin  < 5 
        Ncheb = 1;
    else
        Ncheb = Nstart;
    end
    minc = chebc(f, leftb, rightb, Ncheb);
    while (abs(minc)>OE)
        Ncheb = Ncheb + 1;
        c = chebc(f, leftb, rightb, Ncheb);
        minc = c(Ncheb);
    end
end