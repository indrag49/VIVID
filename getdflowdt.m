function rho = getdflowdt(t,x0,y0,t0,parVec)
b = parVec(1);
A = [0, 1; -1, -b];
rho = A*expm((t-t0)*A)*([x0;y0] - getPhi(t0,parVec)) + getdPhidt(t,parVec);
