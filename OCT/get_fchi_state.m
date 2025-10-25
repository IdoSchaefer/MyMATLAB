function ftarget = get_fchi_state(target)
    ftarget = @(psiT) (target'*psiT)*target;
end