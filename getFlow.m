function rho = getFlow(t,x0,y0,t0,parVec)
b = parVec(1);
omega = parVec(2);
F = parVec(3);
M = [0, 1; -1, -b];
rho = expm((t-t0)*M)*([x0;y0] - S(t0,b,omega,F)) + S(t,b,omega,F);
