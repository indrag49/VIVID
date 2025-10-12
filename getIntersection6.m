function [t3,rho] = getIntersection6(x0,y0,t0,parVec,tEvol,tOff)
iterMax = 40;
tol = 1e-10;
parVec2 = [parVec,x0,y0,t0,2];
t3 = bisectGeneral('flowComponent',t0+tEvol-tOff,t0+tEvol+tOff,iterMax,tol,parVec2);
rho = getFlow(t3,x0,y0,t0,parVec);
