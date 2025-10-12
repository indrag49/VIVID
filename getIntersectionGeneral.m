function [t3,rho] = getIntersectionGeneral(x0,y0,t0,parVec,dt,tDir,runType,xThreshold)

iterMax = 40;
tol = 1e-10;

tOld = t0;
hasCrossedMan = 0;
if runType == 1
    compValue = 1;
elseif runType == 2
    compValue = 2;
elseif runType == 3
    compValue = 2;
end
yValNew = y0;
while hasCrossedMan == 0
    tNew = tOld + tDir*dt;
    rho = getFlow(tNew,x0,y0,t0,parVec);
    if runType == 1
        % manifold is x=0
        myBoolean = rho(1) < 0;
    elseif runType == 2
        % manifold is y=0
        yValOld = yValNew;
        yValNew = rho(2);
        myBoolean = yValOld*yValNew < 0;
    elseif runType == 3
        % manifold is y=0, with x > xThreshold
        yValOld = yValNew;
        xValNew = rho(1);
        yValNew = rho(2);
        myBoolean = (yValOld > 0 && yValNew < 0 && xValNew > xThreshold);
    end
    if myBoolean
        hasCrossedMan = 1;
    else
        tOld = tNew;
    end
end
parVec2 = [parVec,x0,y0,t0,compValue];
t3 = bisectGeneral('flowComponent',tOld,tNew,iterMax,tol,parVec2);
rho = getFlow(t3,x0,y0,t0,parVec);
