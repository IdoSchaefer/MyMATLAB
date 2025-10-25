function result = Pnodes(x, nodes)
% The function computes a polynomial with zeros specified by the variable
% nodes. The leading coefficient is 1.
% x: the values in which the polynomial is to be evaluated.
    Nnodes = length(nodes);
    sz = size(x);
    result = ones(sz);
    for inodes = 1:Nnodes
        result = result.*(x - nodes(inodes));
    end
end