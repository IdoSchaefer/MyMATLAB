function E=plotfrob(a,b)
%E is the 1 minus the normalised Frobenius scalar product between matrices
%  A = a(i) and B = b(i) for i running on all time points.
ll=size(a);ll=ll(3);
E=zeros(1,ll);

for i=1:ll
    A = a(:,:,i);
    B = b(:,:,i);
    E(i) = frob(A,B)/sqrt(frob(A,A)*frob(B,B));
end
E=1-E;

end