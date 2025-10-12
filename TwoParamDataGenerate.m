function [param2Vec,data_all] = TwoParamDataGenerate

runNum = 22;

if runNum == 21
    pp = 2;
    b = 0.04;
    r = .9;
    param2Type = 2;         % omega
    param2 = 0.81;
    param2OtherSteps = -.0001;
    param2FirstStep = param2OtherSteps/10;
    yOff = .001;
    contNumMax = 1200;
    critEig = -1;           %-1 -- PD;
    yGuess = 0.065296175795362;
    FGuess = 0.347533142901875;
    tGuess = 7.853623325056863;

    omega = param2;
    tEvolGuess = pp*2*pi/omega;
elseif runNum == 22
    pp = 2;
    b = 0.04;
    r = .9;
    param2Type = 2;         % omega
    param2 = 0.81;
    param2OtherSteps = .0001;
    param2FirstStep = param2OtherSteps/10;
    yOff = .001;
    contNumMax = 1900;
    critEig = -1;           %-1 -- PD;
    yGuess = 0.065296175795362;
    FGuess = 0.347533142901875;
    tGuess = 7.853623325056863;

    omega = param2;
    tEvolGuess = pp*2*pi/omega;
elseif runNum == 23
    pp = 2;
    b = 0.04;
    r = .9;
    param2Type = 2;         % omega
    param2 = 0.799;
    param2OtherSteps = .00001;
    param2FirstStep = param2OtherSteps/10;
    yOff = .01;
    contNumMax = 5000;
    critEig = 1;           %1 -- SN;
    yGuess = 0.124145165087902;
    FGuess = 0.362836291689990;
    tGuess = 0.071735875674158;

    omega = param2;
    tEvolGuess = pp*2*pi/omega;
elseif runNum == 24
    pp = 2;
    b = 0.04;
    r = .9;
    param2Type = 2;         % omega
    param2 = 0.799;
    param2OtherSteps = -.00002;
    param2FirstStep = param2OtherSteps/10;
    yOff = .003;
    contNumMax = 545;
    critEig = 1;           %1 -- SN;
    yGuess = 0.124145165087902;
    FGuess = 0.362836291689990;
    tGuess = 0.071735875674158;
    omega = param2;
    tEvolGuess = pp*2*pi/omega;
end


iterMaxBisection = 40;
tolBisection = 1e-5;

NewtonOk = 1;
contNum = 1;

param2Vec = zeros(contNumMax,1);
data_all = zeros(contNumMax,12); % [F,e1,e2,x0,t0,y1,t1,y2,t2,x3,t3,tEvol]
F_grazVec = zeros(contNumMax,1);

tic
while contNum <= contNumMax && NewtonOk == 1
    if param2Type == 2
        omega = param2;
    end
    
    if contNum >= 3
        yGuess = 2*data_all(contNum-1,6) - data_all(contNum-2,6);
        FGuess = 2*data_all(contNum-1,1) - data_all(contNum-2,1);
        tGuess = 2*data_all(contNum-1,7) - data_all(contNum-2,7);
        tEvolGuess = 2*data_all(contNum-1,12) - data_all(contNum-2,12);
    end
    y_a = yGuess - yOff;
    y_b = yGuess + yOff;
    
    parVec3 = [FGuess,tGuess,tEvolGuess,b,omega,r,critEig];
    [~,NewtonData,flg] = bisectGeneralNew('fB6',y_a,y_b,iterMaxBisection,tolBisection,parVec3);

    if flg == 0
        NewtonOk = 0;
        param2Vec = param2Vec(1:contNum-1,:);
        data_all = data_all(1:contNum-1,:);
        F_grazVec = F_grazVec(1:contNum-1,:);
    else
        param2Vec(contNum) = param2;
        data_all(contNum,:) = NewtonData;
        F_grazVec(contNum) = sqrt((1-omega^2)^2 + (b*omega)^2);
        disp([num2str(contNum),' / ',num2str(contNumMax)])
        if contNum == 1
            param2 = param2 + param2FirstStep;
        else
            param2 = param2 + param2OtherSteps;
        end
        contNum = contNum + 1;
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(67)
clf
hold on

plot(data_all(:,1),param2Vec,'bx-', 'LineWidth', 1)
plot(F_grazVec,param2Vec,'gx-', 'LineWidth', 1)

DAT = [data_all, param2Vec];

if runNum==21
    dlmwrite('pp2_omega0.8_PD_testing.csv', DAT, 'delimiter', ',', 'precision', '%.16f')
elseif runNum==22
    dlmwrite('pp2_omega0.8_PDb_testing.csv', DAT, 'delimiter', ',', 'precision', '%.16f')
elseif runNum==23
    dlmwrite('pp2_omega0.8_SN_testing.csv', DAT, 'delimiter', ',', 'precision', '%.16f')
elseif runNum==24
    dlmwrite('pp2_omega0.8_SNb_testing.csv', DAT, 'delimiter', ',', 'precision', '%.16f')
end