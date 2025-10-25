function meanabsp = meanproj(psi, wavefcns)
% The function computes the mean projections in time, of a number of wave
% functions, on a series of wave function values in different times.
% input:
% psi: contains the wave function values in all the time points, in
% different columns.
% wavefcns: a matrix. Contains the wave functions to project onto, in
% differnt columns.
% output:
% meanabsp: a row vector. The i'th term is the maximal absolute value of the
% projections of the wave function on the i'th column of wavefcns on psi in
% all times.
    Nwf = size(wavefcns, 2);
    meanabsp = zeros(1, Nwf);
    for wfi = 1:Nwf
        alltproj = wavefcns(:, wfi)'*psi;
        meanabsp(wfi) = mean(abs(alltproj));
    end
end