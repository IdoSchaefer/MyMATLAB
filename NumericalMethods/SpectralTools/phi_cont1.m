function phic = phi_cont1(phi)
% The function gets a non-continuous instantaneous phase vector phi, 
% and returns a continuous phase vector phic.
% The function handles with discontinuities of integer multiplies of pi/2.
    Nphi = length(phi);
    dphi = phi(2:Nphi) - phi(1:(Nphi - 1));
    pi_half_jumps = (abs(dphi)>0.45*pi).*round(dphi/(pi/2));
    fixjumps = 0;
    phic = phi;        
    for iphi = 2:Nphi
        fixjumps = fixjumps - pi_half_jumps(iphi - 1)*pi/2;
        phic(iphi) = phic(iphi) + fixjumps;
    end
end