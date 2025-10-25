function [pout, pouti] = Lejap1(oldp, newp, Nout, capacity)
% The function returns the Nout Leja points, chosen from the set of points:
% newp. The old point set is oldp.
% capacity: The capacity of the approximation domain; required for stability. 
    if nargin<4
        capacity = 1;
    end
    Nold = length(oldp);
    dim = size(newp);
    if dim(1) == 1
        Nnew = dim(2);
        pout = zeros(1, Nout);
        difmul = ones(1, Nnew);
        if nargout>1
            pouti = zeros(1, Nout);
        end
    else
        Nnew = dim(1);
        pout = zeros(Nout, 1);
        difmul = ones(Nnew, 1);
        if nargout>1
            pouti = zeros(Nout, 1);
        end
    end
    is_unused = true(Nnew, 1);
    outi = 1;
% The case of an empty oldp requires a special treatment for the 1'st
% point - we have to choose it smartly:
    if Nold == 0        
        [~, maxi] = max(abs(newp));
        pout(outi) = newp(maxi);
        is_unused(maxi) = false;
        difmul(maxi) = 0;
        if nargout>1
            pouti(outi) = maxi;
        end
        Nold = 1;
        oldp = pout(outi);
        outi = outi + 1;
    end
    for oldi = 1:Nold
        difmul(is_unused) = difmul(is_unused).*abs(newp(is_unused) - oldp(oldi))/capacity;
    end
    for outi = outi:Nout
        [~, maxi] = max(difmul);
        pout(outi) = newp(maxi);
        is_unused(maxi) = false;
        difmul(maxi) = 0;
        if nargin>1
            pouti(outi) = maxi;
        end
        difmul(is_unused) = difmul(is_unused).*abs(newp(is_unused) - pout(outi))/capacity;
    end
end    