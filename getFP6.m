function NewtonData = getFP6(y,parVec3)

tolNewton = 1e-10;
maxNewtonSteps = 10;
stepNum = 1;

[res,DGdyt,dGdF,Dfp,fullData,flg] = evaluateG6(y,parVec3);

if flg == 1
    DGdFt = [dGdF,DGdyt(:,2)];
    e = sort(eig(Dfp));
    while flg == 1 && stepNum <= maxNewtonSteps && norm(res) > tolNewton
        FtNew = [parVec3(1);parVec3(2)] - DGdFt\res;

        parVec3(1) = FtNew(1);
        parVec3(2) = FtNew(2);
        [res,DGdyt,dGdF,Dfp,fullData,flg] = evaluateG6(y,parVec3);
        if flg == 1
            DGdFt = [dGdF,DGdyt(:,2)];
            e = sort(eig(Dfp));
            stepNum = stepNum + 1;
        end
    end
end

if flg == 1
    NewtonData = [parVec3(1),e(1),e(2),fullData];
else
    NewtonData = [];
end

