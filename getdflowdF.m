function rho = getdflowdF(t,x0,y0,t0,parVec)
b = parVec(1);
A = [0, 1; -1, -b];
rho = getdPhidF(t,parVec) - expm((t-t0)*A)*getdPhidF(t0,parVec);
