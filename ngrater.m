function vi = ngrater(v, x)
%     vi = 1;
%     while v(vi) < x
%         vi = vi + 1;
%     end
    N = length(v);
    bigi = N;
    smalli = 1;
    vi = floor((bigi + smalli)/2);
    while v(vi)~=x && (bigi-smalli)>1
        if v(vi) > x
            bigi = vi;
%         elseif v(vi) < x
        else
            smalli = vi;
%         else
%             break
        end
        vi = floor((bigi + smalli)/2);
    end
    if v(vi)<x
        vi = vi + 1;
    end
end