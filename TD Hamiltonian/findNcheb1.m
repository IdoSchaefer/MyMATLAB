function Ncheb = findNcheb1(f, leftb, rightb, relE, x)
% The function finds the minimal number of terms (Ncheb) in a Chebyshev expansion
% of a function f, in the region [leftb, rightb], needed to get a relative
% error below relE, in the point x. Without an input for x, the function
% takes a default value.
    if nargin<5
        if abs(f(leftb))> abs(f(rightb));
            x = leftb;
        else
            x = rightb;
        end
    end
    Ncheb = 2;
    while(abs(Echeb(f, x, leftb, rightb, Ncheb))>relE && Ncheb<1000)
        Ncheb = Ncheb + 1;
    end
end
