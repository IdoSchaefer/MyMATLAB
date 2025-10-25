function [vnew, origi] = rand_order(v)
% The function returns a vector vnew, in which the order of the terms in
% the vector v is randomized.
% origi: The original indices of the terms of vnew in v.
%%%% vi: The new indices of the terms of v in vnew. vnew(vi) = v.
%%%% I'm not sure that this will be useful, so I commented the relevant lines.
    N = length(v);
    [~, origi] = sort(rand(1, N));
    vnew = v(origi);
%    vi(newi) = 1:N;
end