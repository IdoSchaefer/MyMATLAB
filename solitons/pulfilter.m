numfil = input('How many pulses do you want to save?');
highpeaks = zeros(5, numfil, toptimei);
nsavedpul = zeros(2, numfil);
for ti = 1:toptimei
    ezer = allpeaks(1:2, :, ti);
    for savedpuli = 1:numfil
        [stam, nhighest] = min(ezer(1, :));
        nsavedpul(1, savedpuli) = nhighest;
        nsavedpul(2, savedpuli) = ezer(2, nhighest);
        ezer(1, nhighest) = r0;
    end
% to sort the pulses according to the position, for easier arrangement,
% with arrpul.
    for i = 1:numfil-1
        for j = 1:numfil-i
            if nsavedpul(2, j+1)<nsavedpul(2, j)
                temp = nsavedpul(:,j);
                nsavedpul(:, j) = nsavedpul(:, j+1);
                nsavedpul(:, j+1) = temp;
            end
        end
    end
    highpeaks(:, : ,ti) = allpeaks(:, nsavedpul(1, :), ti);
end