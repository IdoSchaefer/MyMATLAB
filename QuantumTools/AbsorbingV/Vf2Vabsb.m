function V = Vf2Vabsb(Vf, Vabs, x)
% The function returns a potential vector V of the form of the function 
% handle @(x) Vf(x), with absorbing boundaries specified by the vector
% Vabs.
    Nx = length(x);
    Nxabs = length(Vabs);
    V = [Vabs(Nxabs:-1:1) + Vf(x(Nabs)); Vf(x((Nabs + 1):(Nx - Nabs))); Vabs + Vf(x(Nx - Nabs + 1))];
end