function z = flowComponent(t,parVec2)
parVec = parVec2(1:3);
x0 = parVec2(4);
y0 = parVec2(5);
t0 = parVec2(6);
compValue = parVec2(7);
rho = getFlow(t,x0,y0,t0,parVec);
z = rho(compValue);
