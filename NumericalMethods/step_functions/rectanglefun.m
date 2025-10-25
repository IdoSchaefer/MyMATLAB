function result = rectanglefun(x, a, b)
% The function computes the rectangle function:
% f(x) = 1, a<=x<=b
% f(x) = 0, otherwise
    result = double(x>=a & x<=b);
%     if x>=a && x<=b
%         result = 1;
%     else
%         result = 0;
%     end
end