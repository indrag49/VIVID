function [res,DGdyt,dGdF,Dfp,fullData,flg] = evaluateG6(y1,parVec3)

F = parVec3(1);
t1 = parVec3(2);
tEvol = parVec3(3);
b = parVec3(4);
omega = parVec3(5);
r = parVec3(6);

xThreshold = -.5;
dt = 2*pi/omega*.05;
tOff = 2*pi/omega*.2; %05;

res = [];
DGdyt = [];
dGdF = [];
Dfp = [];
fullData = [];

timeDir = sign(y1);
A = [0,1;-1,-b];
parVec = [b,omega,F];
flg = 1;

% compute x0, t0 and derivative
[t0,~] = getIntersectionGeneral(0,y1,t1,parVec,dt,timeDir,2,xThreshold); 
if flg == 0
    disp('ERROR: x0, t0, computation failed')
else
    rho = getFlow(t0,0,y1,t1,parVec);
    x0 = rho(1);
    Dflow = expm((t0-t1)*A);
    uuu = getdflowdt(t0,0,y1,t1,parVec);
    vvv = getdflowdt0(t0,0,y1,t1,parVec);
    www = getdflowdF(t0,0,y1,t1,parVec);
    M0 = 1/uuu(2)*[Dflow(1,2)*uuu(2)-Dflow(2,2)*uuu(1), uuu(2)*vvv(1)-uuu(1)*vvv(2);-Dflow(2,2),-vvv(2)];
    M1 = inv(M0);
    N0 = [-uuu(1)*www(2)/uuu(2)+www(1);-www(2)/uuu(2)];

    % compute y2, t2 and derivative
    y2 = -r*y1;
    t2 = t1;
    M2 = [-r,0;0,1];

    % compute x3, t3 and derivative
    [t3,~] = getIntersectionGeneral(0,y2,t2,parVec,dt,-timeDir,2,xThreshold);
    if flg == 0
        disp('ERROR: x3, t3, computation failed')
    else
        rho = getFlow(t3,0,y2,t2,parVec);
        x3 = rho(1);
        Dflow = expm((t3-t2)*A);
        uuu = getdflowdt(t3,0,y2,t2,parVec);
        vvv = getdflowdt0(t3,0,y2,t2,parVec);
        www = getdflowdF(t3,0,y2,t2,parVec);
        M3 = 1/uuu(2)*[Dflow(1,2)*uuu(2)-Dflow(2,2)*uuu(1), uuu(2)*vvv(1)-uuu(1)*vvv(2); -Dflow(2,2), -vvv(2)];
        N3 = [-uuu(1)*www(2)/uuu(2)+www(1);-www(2)/uuu(2)];

        % computes p^th iterate of the global map and derivative
        [t4,~] = getIntersection6(x3,0,t3,parVec,tEvol,tOff);
        if flg == 0
            disp('ERROR: x4, t4, computation failed')
        else
            Dflow = expm((t4-t3)*A);
            uuu = getdflowdt(t4,x3,0,t3,parVec);
            vvv = getdflowdt0(t4,x3,0,t3,parVec);
            www = getdflowdF(t4,x3,0,t3,parVec);
            M4 = 1/uuu(2)*[Dflow(1,1)*uuu(2) - Dflow(2,1)*uuu(1), uuu(2)*vvv(1) - uuu(1)*vvv(2); -Dflow(2,1), -vvv(2)];
            N4 = [-uuu(1)*www(2)/uuu(2)+www(1);-www(2)/uuu(2)];
            rho = getFlow(t4,x3,0,t3,parVec);
            x4 = rho(1);
            fullData = [x0,t0,y1,t1,y2,t2,x3,t3,t4-t3];
            t4 = mod(t4,2*pi/omega);
        end

        if flg == 1
            tDiff = clockDifference(t0,t4,2*pi/omega);
            res = [x4-x0;tDiff];
            DGdyt = M4*M3*M2 - M0;
            dGdF = M4*N3 + N4 - N0;
            Dfp = M4*M3*M2*M1;

            
        end
    end
end
