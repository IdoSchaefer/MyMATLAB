function phic = phi_cont(phi)
% The function gets a non-continuous instantaneous phase vector phi, 
% and returns a continuous phase vector phic.
    Nphi = length(phi);
    dphi = phi(2:Nphi) - phi(1:(Nphi - 1));
    % Finding the indices of the positive jumps in the phase:
    pos_jumpi = find(dphi>pi);
    Npos = length(pos_jumpi);
    % Finding the indices of the negative jumps in the phase:
    neg_jumpi = find(dphi<-pi);
    Nneg = length(neg_jumpi);
    phic = phi;
    % Fixing the positive jumps:
    for posi = 1:Npos
        phic((pos_jumpi(posi) + 1):Nphi) = phic((pos_jumpi(posi) + 1):Nphi) - 2*pi;
    end
    % Fixing the negative jumps:
    for negi = 1:Nneg
        phic((neg_jumpi(negi) + 1):Nphi) = phic((neg_jumpi(negi) + 1):Nphi) + 2*pi;
    end
end