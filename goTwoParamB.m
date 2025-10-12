function goTwoParamB

figure(83)
set(gcf,'Position',[100 100 500 440])
set(gca,'Position',[0.14,0.14,0.8,0.8])
clf
hold on

xMin = 0.34;
xMax = 0.37;
yMin = 0.785;
yMax = 0.82;

xlim([xMin xMax]);
ylim([yMin yMax]);

b = 0.04;
r = .9;

lw=5;
sz=18;

OO_ori = linspace(yMin, yMax, 10000);
FF_ori = sqrt((1 - OO_ori.^2).^2 + b^2 .* OO_ori.^2);

plot(FF_ori, OO_ori, 'g-', 'LineWidth', lw)

% Plot grazing curves

h = 1e-4;
Gtol = 1e-5;
numNewtonMax = 10;

omegaStart = 0.839;
omegaStep = .002;
numStep =81;
FtGuess = [.5214;1.];
impactVec = [0,1];
xThreshold = -.5;
dtFrac = .05;

omegaVec = zeros(numStep,1);
FtAll = zeros(numStep,2);
for k = 1:numStep
    if k == 1
        omega = omegaStart;
    else
        omega = omega + omegaStep;
    end
    if k >= 3
        FtGuess = 2*FtAll(k-1,:)' - FtAll(k-2,:)';
    end
    dt = 2*pi/omega*dtFrac;
    parVec7 = [b,omega,r,xThreshold,dt,0];

    % Newton!!!
    j = 1;
    G = getAuxGrazFunc(FtGuess(1),FtGuess(2),impactVec,parVec7);
    while norm(G) > Gtol && j <= numNewtonMax
        Gpa = getAuxGrazFunc(FtGuess(1)+h,FtGuess(2),impactVec,parVec7);
        Gma = getAuxGrazFunc(FtGuess(1)-h,FtGuess(2),impactVec,parVec7);
        Gpb = getAuxGrazFunc(FtGuess(1),FtGuess(2)+h,impactVec,parVec7);
        Gmb = getAuxGrazFunc(FtGuess(1),FtGuess(2)-h,impactVec,parVec7);
        DG = 1/(2*h)*[Gpa-Gma,Gpb-Gmb];
        FtGuess = FtGuess - DG\G;
        G = getAuxGrazFunc(FtGuess(1),FtGuess(2),impactVec,parVec7);
        j = j+1;
    end

    omegaVec(k) = omega;
    FtAll(k,:) = FtGuess';
    disp(k)
end
plot(FtAll(:,1), omegaVec, 'Color', '#800080', 'LineWidth', lw); % Dark purple

Purp1A = FtAll(:,1);
Purp1B = omegaVec;


h = 1e-4;
Gtol = 1e-5;
numNewtonMax = 10;

omegaStart = 0.839;
omegaStep = -.0002;
numStep =255;
FtGuess = [.5214;1.];
impactVec = [0,1];
xThreshold = -.5;
dtFrac = .05;

omegaVec = zeros(numStep,1);
FtAll = zeros(numStep,2);
for k = 1:numStep
    if k == 1
        omega = omegaStart;
    else
        omega = omega + omegaStep;
    end
    if k >= 3
        FtGuess = 2*FtAll(k-1,:)' - FtAll(k-2,:)';
    end
    dt = 2*pi/omega*dtFrac;
    parVec7 = [b,omega,r,xThreshold,dt,0];

    % Newton!!!
    j = 1;
    G = getAuxGrazFunc(FtGuess(1),FtGuess(2),impactVec,parVec7);
    while norm(G) > Gtol && j <= numNewtonMax
        Gpa = getAuxGrazFunc(FtGuess(1)+h,FtGuess(2),impactVec,parVec7);
        Gma = getAuxGrazFunc(FtGuess(1)-h,FtGuess(2),impactVec,parVec7);
        Gpb = getAuxGrazFunc(FtGuess(1),FtGuess(2)+h,impactVec,parVec7);
        Gmb = getAuxGrazFunc(FtGuess(1),FtGuess(2)-h,impactVec,parVec7);
        DG = 1/(2*h)*[Gpa-Gma,Gpb-Gmb];
        FtGuess = FtGuess - DG\G;
        G = getAuxGrazFunc(FtGuess(1),FtGuess(2),impactVec,parVec7);
        j = j+1;
    end

    omegaVec(k) = omega;
    FtAll(k,:) = FtGuess';
    disp(k)
end

plot(FtAll(:,1), omegaVec, 'Color', '#800080', 'LineWidth', lw);

Purp2A = FtAll(:,1);
Purp2B = omegaVec;

xi = sqrt(1 - b^2/4);
phi = r;

% for omega \approx 0.8
p = 2;
omega = xi/(5/4) %p=2
F_graz = sqrt(omega^4+(b^2-2)*omega^2+1)

plot([F_graz 0.358], [omega 0.7985], 'k-', 'LineWidth', 1)
text(0.3525, 0.7979, 'resonance', 'FontSize', 20, 'Interpreter', 'latex')
plot(sqrt((1 - omega^2).^2 + b^2*omega^2), omega, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 11)

DATA_PD = readmatrix('pp2_omega0.8_PD.csv');
FVec_PD_data = DATA_PD(:,1);
omegaVec_PD_data = DATA_PD(:, end);
plot(FVec_PD_data,omegaVec_PD_data, '-', 'color', 'black', 'LineWidth', lw)
plot([FVec_PD_data(end) 0.361674],[omegaVec_PD_data(end) 0.79984], '-', 'color', 'black', 'LineWidth', lw)

% xFill = [FVec_PD_data(1:end-35); flipud(Purp2A(1:end-35))];
% yFill = [omegaVec_PD_data(1:end-35); flipud(Purp2B(1:end-35))];
% fill(xFill, yFill, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

DATA_PD2 = readmatrix('pp2_omega0.8_PDb.csv');
FVec_PD2_data = DATA_PD2(:,1);
omegaVec_PD2_data = DATA_PD2(:, end);
plot(FVec_PD2_data,omegaVec_PD2_data, '-', 'color', 'black', 'LineWidth', lw)


% xFill = [FVec_PD2_data; flipud(Purp1A)];
% yFill = [omegaVec_PD2_data; flipud(Purp1B)];
% fill(xFill, yFill, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

DATA_SN = readmatrix('pp2_omega0.8_SN.csv');
FVec_SN_data = DATA_SN(:,1);
omegaVec_SN_data = DATA_SN(:, end);
plot(FVec_SN_data,omegaVec_SN_data, '-', 'color', 'black', 'LineWidth', lw)

DATA_SN2 = readmatrix('pp2_omega0.8_SNb.csv');
FVec_SN2_data = DATA_SN2(:,1);
omegaVec_SN2_data = DATA_SN2(:, end);
plot(FVec_SN2_data,omegaVec_SN2_data, '-', 'color', 'black', 'LineWidth', lw)

% disp(FVec_SN2_data(end))
% disp(omegaVec_SN2_data(end))

xFill = [Purp2A(end:-1:1);Purp1A;FVec_PD2_data(end:-1:1);FVec_PD_data;FVec_SN_data(end:-1:1);FVec_SN2_data];
yFill = [Purp2B(end:-1:1);Purp1B;omegaVec_PD2_data(end:-1:1);omegaVec_PD_data;omegaVec_SN_data(end:-1:1);omegaVec_SN2_data];
fill(xFill, yFill, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('$\mathcal{A}$', 'Interpreter', 'latex', 'FontSize',30);
ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 30, 'Rotation', 0);

% yticks([0.78, 0.81, 0.84, 0.86]);
 
xticks([0.34, 0.35, 0.36, 0.37]);
xticklabels({'0.34', '0.35', '0.36', '0.37'});  
yticks([0.79, 0.8, 0.81, 0.82]);
yticklabels({'0.79', '0.8', '0.81', '0.82'});  

set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);

text(0.364, 0.7995, '$\mathcal{A}_{\rm graz}$', 'Interpreter', 'latex', 'FontSize',30)
text(0.361, 0.7905, 'secondary', 'FontSize', 20, 'Interpreter', 'latex');
text(0.361, 0.7888, 'grazing', 'FontSize', 20, 'Interpreter', 'latex');
text(0.3458, 0.8125, 'PD', 'Interpreter', 'latex', 'FontSize',20)
text(0.3617, 0.796, 'SN', 'Interpreter', 'latex', 'FontSize',20)

plot(0.350751, 0.788158, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 11)

function G = getAuxGrazFunc(F,t0,impactVec,parVec7)

b = parVec7(1);
omega = parVec7(2);
r = parVec7(3);
xThreshold = parVec7(4);
dt = parVec7(5);
plotStuff = parVec7(6);

parVec = [b,omega,F];

if plotStuff == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(102)
    clf
    hold on
    grid on
end

tStart = t0;
x0 = 0;
for oscNum = 1:length(impactVec)
    if impactVec(oscNum) == 1
        % evaluate disc map
        [t1,rho1] = getIntersectionGeneral(x0,0,t0,parVec,dt,-1,1,xThreshold);
        t2 = t1;
        y2 = -r*rho1(2);
        [t3,rho3] = getIntersectionGeneral(0,y2,t2,parVec,dt,-1,2,xThreshold);
        x3 = rho3(1);
    else
        t3 = t0;
        x3 = x0;
    end
    % evaluate global map
    [t4,rho4] = getIntersectionGeneral(x3,0,t3,parVec,dt,1,3,xThreshold);

    if plotStuff == 1
        ttt = linspace(t3,t4,500);
        zzz = zeros(2,length(ttt));
        for j = 1:length(ttt)
            zzz(:,j) = getFlow(ttt(j),x3,0,t3,parVec);
        end
        plot(zzz(1,:),zzz(2,:),'kx-')
    end

    t0 = t4;
    x0 = rho4(1);
end
G = [x0;clockDifference(tStart,t0,2*pi/omega)];
