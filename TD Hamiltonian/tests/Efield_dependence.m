load qubits2HOtest_problem Hu2modesa H1u2modesa
Efield1 = zeros(19, 101);
for fieldi=1:101
    Efield1(:,fieldi) = eig(Hu2modesa - (fieldi-1)*0.01*H1u2modesa);
end
