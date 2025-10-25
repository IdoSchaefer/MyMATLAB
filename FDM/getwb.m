U0b = U0(jj:Kb + jj - 1, jj:Kb + jj - 1);
U1b = U1(jj:Kb + jj - 1, jj:Kb + jj - 1);
[Bb, Dub] = eig(U1b, U0b);
wb = -log(diag(Dub).')/(i*dt);
wb(real(wb)<0) = wb(real(wb)<0) + 2*pi/dt;