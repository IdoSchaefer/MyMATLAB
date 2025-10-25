function result = sqnorm(psi)
% The function returns the square of the norm of psi, <psi|psi>.
% Suitable also for psi in many time steps, each represented in a column.
    result = sum(conj(psi).*psi);
end