function dphidt = getdPhidt(tt,parVec)
b = parVec(1);
omega = parVec(2);
F = parVec(3);
dphidt = omega*F/(omega^4+(b^2-2)*omega^2+1)*[b*omega,-(1-omega^2);-omega*(1-omega^2),-b*omega^2]*[cos(omega*tt);sin(omega*tt)];
