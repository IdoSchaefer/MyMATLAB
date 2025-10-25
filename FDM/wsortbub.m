for i1 = 1:Kb - 1
    for i2 = 1:(Kb - i1)
        if real(wb(i2 + 1))<real(wb(i2))
            ezer = wb(i2);
            wb(i2) = wb(i2 + 1);
            wb(i2 + 1) = ezer;
            ezerB = Bb(:, i2);
            Bb(:, i2) = Bb(:, i2 + 1);
            Bb(:, i2 + 1) = ezerB;
        end
    end
end
