A = [5 4 3 2 1];
for i = 1:4
    for j = 1:(5 - i)
        if A(j+1)<A(j)
            ezer = A(j);
            A(j) = A(j+1);
            A(j+1) = ezer;
        end
    end
end
A