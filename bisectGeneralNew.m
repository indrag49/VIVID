% general bisection method
function [a3,DATA,flg] = bisectGeneralNew(func,a1,a2,iterMax,tol,parVec)

[z1,~,flg1] = feval(func,a1,parVec);
[z2,~,flg2] = feval(func,a2,parVec);
if flg1 == 0 || flg2 == 0
    a3 = [];
    DATA = [];
    flg = 0;
else
    if sign(z1) == sign(z2)
        flg = 0;
        disp('ERROR: initial interval invalid for bisection method')
        disp(['a1: ',num2str(a1),'  z1: ',num2str(z1)])
        disp(['a2: ',num2str(a2),'  z2: ',num2str(z2)])
        a3 = [];
        DATA = [];
    else
        a3 = (a1+a2)/2;
        [z3,DATA,flg3] = feval(func,a3,parVec);
        j = 1;
        while flg3 == 1 && j <= iterMax && abs(z3) > tol
            if sign(z3) == sign(z1)
                a1 = a3;
                z1 = z3;
            else
                a2 = a3;
            end
            a3 = (a1+a2)/2;
            [z3,DATA,flg3] = feval(func,a3,parVec);
            j = j+1;
        end
        if j > iterMax
            disp('WARNING: iterMax exceeded')
        end
        if flg3 == 0
            a3 = [];
            DATA = [];
            flg = 0;
        else            
            flg = 1;
        end
    end
end
