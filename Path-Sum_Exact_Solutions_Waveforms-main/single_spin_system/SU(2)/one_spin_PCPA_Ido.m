function [rho,elapsedTime] = one_spin_PCPA_Ido(para,TMAX,NT)
%% M. Foroozandeh, P.-L. Giscard, 04/2022
% para : parameters as set by paragen
% TMAX : in second, maximum evolution time
% NT : number of evaluation points
% Time evolution of the elements of density matrix for a single spin-1/2
% propagation by PCPA (piecewise-constant propagator approximation) expansion

t_array = linspace(0,TMAX,NT);
tres = t_array(2);

smfactor = para.n;
offs_f = para.deltaf;
bandwidth = para.DeltaF;
phi0 = para.Phi0;
tau_p = para.taup;
Omega = para.Omega;
omega1 = para.omega1;
offs_t = para.deltat;

tic;

%pauli matrices

Sigmax = 0.5*[0,1;1,0];Sigmay = 0.5*[0,-1i;1i,0];Sigmaz = 0.5*[1,0;0,-1];

% building multi-state operators for both spins
Lx = Sigmax;
Ly = Sigmay;
Lz = Sigmaz;

Omega_L = Omega;

% pulse

waveform=chirp_fun(tau_p,bandwidth,phi0,omega1,offs_t,offs_f,smfactor,t_array);

rho_0 = 2*Lz; % initial state

% Preallocating memory (Ido):

rho = zeros(2, 2, length(waveform));

% This takes the offset and pulse information and run the numerical
% simulation and then plot the time evolution of the elements of the final density matrix

for i=1:length(waveform)
    
    H = Omega_L*Lz + real(waveform(i))*Lx + imag(waveform(i))*Ly;
    
    U = expm(-1i*tres*H);
    
    %rho(:,:,i) =expm(-1i*tres*H)*rho_0*expm(1i*tres*H); %  (eye(2)+(-1i)*tres*H{i})*rho_0*(eye(2)+(1i)*tres*H{i});%
    rho(:,:,i) =U*rho_0*U';
    rho_0=rho(:,:,i);
    
    
end

elapsedTime = toc;


end

function pulse=chirp_fun(tau_p,bandwidth,phi0,omega1,offs_t,offs_f,smfactor,t_array)

Cx = (exp(-(2^(smfactor+2))*((t_array-offs_t)/tau_p).^smfactor)).*(omega1*cos(phi0+(pi*bandwidth*((t_array-offs_t).^2)/tau_p)-2*pi*offs_f*(t_array-offs_t)));
Cy = (exp(-(2^(smfactor+2))*((t_array-offs_t)/tau_p).^smfactor)).*(omega1*sin(phi0+(pi*bandwidth*((t_array-offs_t).^2)/tau_p)-2*pi*offs_f*(t_array-offs_t)));

pulse = complex(Cx,Cy);

end
