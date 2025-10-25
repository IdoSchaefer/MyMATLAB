function ftarget = phi2ftarget(phi)
    ftarget = @(psi) phi*(phi'*psi);
end