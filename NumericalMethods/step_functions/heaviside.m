function result = heaviside(x)
% The Heaviside step function:
    result = double(x>=0);
%     if x<0
%         result = 0;
%     else
%         result = 1;
%     end
end