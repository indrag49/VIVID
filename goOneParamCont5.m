function data_all = goOneParamCont5
%function goOneParamCont5

lw = 4;
ms = 8;

runNum = 21;

if runNum == 21
    pp = 2;
    b = 0.04;
    r = 0.9;
    omega = .81; %PD
    yStart = 0.409097036428900 ;   % impact velocity
    yStep = -.0001;
    FGuessStart = 0.36;
    tGuessStart = 7.755681336950062;   % impact time
    contNumMax = 5000;
    findPD = 1;
        
 elseif runNum == 22
    pp = 2;
    b = 0.04;
    r = 0.9;
    omega = .799; %SN
    yStart = 0.245022580357925 ;   % impact velocity
    yStep = -.0001;
    FGuessStart = 0.363;
    tGuessStart = 0.034154055932653;   % impact time
    contNumMax = 4000;
    findPD = 0;
 elseif runNum == 23
    pp = 2;
    b = 0.04;
    r = 0.9;
    omega = .799; %SN
    yStart = 0.245022580357925 ;   % impact velocity
    yStep = .0001;
    FGuessStart = 0.363;
    tGuessStart = 0.034154055932653;   % impact time
    contNumMax = 2000;
    findPD = 0;
end

F_graz = sqrt(omega^4+(b^2-2)*omega^2+1)



%delta = exp(-2*pi/omega*b)
%tau = 2*exp(-pi/omega*b)*cos(2*pi/omega*sqrt(1-b^2/4))

data_all = zeros(contNumMax,2*pp+11);  % [F,e1,e2,x0,t0,y1,t1,y2,t2,x3,t3,x4i,t4i]

NewtonOk = 1;
contNum = 1;
while contNum <= contNumMax && NewtonOk == 1
    if contNum == 1
        y = yStart;
        FGuess = FGuessStart;
        tGuess = tGuessStart;
    elseif contNum == 2
        FGuess = data_all(1,1);
        tGuess = data_all(1,7);
    else
        FGuess = 2*data_all(contNum-1,1) - data_all(contNum-2,1);
        tGuess = 2*data_all(contNum-1,7) - data_all(contNum-2,7);
    end

    parVec3 = [FGuess,tGuess,pp,b,omega,r];
    NewtonData = getFP5(y,parVec3);
    
    if isempty(NewtonData)
        NewtonOk = 0;
        data_all = data_all(1:contNum-1,:);
    else        
        data_all(contNum,:) = NewtonData;
        disp([num2str(contNum),' / ',num2str(contNumMax)])
        contNum = contNum + 1;
        y = y + yStep;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%



if runNum == 21
    figure(61)
    set(gcf,'Position',[100 100 500 440])
    set(gca,'Position',[0.14,0.14,0.8,0.8])
    clf
    hold on
    grid off

    
    set(groot, 'defaultLegendInterpreter', 'latex');
    
    % annotation(gcf,'arrow',[0.833333333333333 0.866852886405959],...
    % [0.630271878646441 0.533255542590432],'LineWidth', 2);

    
    % Add $\mathcal{A}_{\rm graz}$ label below the x-axis
    text(0.345, -0.8, '$\mathcal{A}_{\rm graz}$', 'FontSize', 25, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
    plot([F_graz 0.3452], [-0.9 -0.85], 'k-', 'LineWidth', 1)

    % text(0.3472, -0.22, 'Complex eigenvalues', 'FontSize', 20, 'Interpreter','latex')
    
    % plot(data_all(:,1),data_all(:,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    % plot(data_all(:,1),data_all(:,3),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_2$') %e2
    % legend('Location', 'southwest');

    % crude approximation to PD from given data
    ind = find(data_all(:,2) < -1,1);

    % plot(data_all(3411:end,1),data_all(3411:end,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    % plot(data_all(3411:end,1),data_all(3411:end,3),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_2$') %e2


    % %filler??
    % plot([0.347627 0.347626], [-0.615201 -0.707883], 'k-', 'LineWidth', lw, 'HandleVisibility', 'off')

    plot(data_all(1:ind,1),data_all(1:ind,6),'b-', 'LineWidth', lw, 'DisplayName', '$y$ (stable)') %y1
    % plot(data_all(ind+1:end,1),data_all(ind+1:end,6),'r-', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    
    indF = find(data_all(:,1) < F_graz,1);
    plot(data_all(ind+1:indF,1),data_all(ind+1:indF,6),'r-', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    plot(data_all(indF+1:end,1),data_all(indF+1:end,6),'r--', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    
    % % plot(data_all(:,1),data_all(:,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    plot(data_all(1:indF,1),data_all(1:indF,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    plot(data_all(indF+1:end,1),data_all(indF+1:end,2),'k--', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    plot(data_all(:,1),data_all(:,3),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_2$') %e2
    
    % plot(data_all(data_all(ind+1:end,1) < F_graz,1), data_all(data_all(ind+1:end,1) < F_graz,6), 'r--', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    % plot(data_all(data_all(ind+1:end,1) >= F_graz,1), data_all(data_all(ind+1:end,1) >= F_graz,6), 'r-', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1

    s = (data_all(ind-1,2)+1)/(data_all(ind-1,2)-data_all(ind,2));
    Fyt_impact_PD = (1-s)*data_all(ind-1,[1,6,7, 13]) + s*data_all(ind,[1,6,7, 13])
    ylim([-1.3, 0.3]);
    xlim([0.344, 0.3482]);
    yline(-1, '-', 'LineWidth', lw-2, 'HandleVisibility', 'off', 'Color', [0.5, 0.5, 0.5])
    xlabel('$\mathcal{A}$', 'Interpreter', 'latex', 'FontSize',30);
    % ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 30, 'Rotation', 0);
    title('$\omega=0.81$', 'FontSize', 30, 'Interpreter','latex');
    plot(Fyt_impact_PD(1),-1,'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [144, 103, 167]/256, 'MarkerSize', ms, 'HandleVisibility', 'off')    
    text(0.34758,-1.06,'PD', 'FontSize', 20)
    % Set selective x-tick values
    xticks([0.344, 0.345, 0.346, 0.347, 0.348]);
    xticklabels({'0.344', '0.345', '0.346', '0.347', '0.348'});  

    text(0.3479, -0.59, '${\rm Re}(\lambda_{1,2})$', 'FontSize', 25, 'Interpreter', 'latex')
    text(0.3445, 0.11, '$\lambda_1$', 'FontSize', 25, 'Interpreter', 'latex')
    text(0.34745, -1.22, '$\lambda_2$', 'FontSize', 25, 'Interpreter', 'latex')
    text(0.346, 0.1, '$y_{\rm imp}$', 'FontSize', 25, 'Interpreter', 'latex')

    xline(F_graz, 'LineStyle', '-', 'Color', 'green', 'LineWidth', lw, 'DisplayName', 'grazing')
    
elseif runNum==22
    figure(62)
    set(gcf,'Position',[100 100 500 440])
    set(gca,'Position',[0.14,0.14,0.8,0.8])
    hold on
    grid off
    
    % set(groot, 'defaultLegendInterpreter', 'latex');
        
    % plot(data_all(:,1),data_all(:,6),'b-', 'LineWidth', lw, 'DisplayName', '$x$') %y1
    % plot(data_all(:,1),data_all(:,7),'r-', 'LineWidth', lw) %t1
    % plot(data_all(:,1),data_all(:,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    % plot(data_all(:,1),data_all(:,3),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_2$') %e2
    % legend('Location', 'southwest');

    %crude approximation to SN from given data
    ind = find(data_all(:,3) > 1,1);

    text(0.36297, -.2, '$\mathcal{A}_{\rm graz}$', 'FontSize', 25, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
    plot([F_graz 0.36299], [-0.32 -0.24], 'k-', 'LineWidth', 1)

    xticks([0.3628, 0.3629, 0.363, 0.3631]);
    xticklabels({'0.3628', '0.3629', '0.363', '0.3631'}); 

    plot(data_all(1:2506,1),data_all(1:2506,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    % plot(data_all(2507:end,1),data_all(2507:end,2),'k--', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    % plot(data_all(1:2506,1),data_all(1:2506,3),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_2$') %e2

    plot(data_all(1:ind,1),data_all(1:ind,6),'b-', 'LineWidth', lw, 'DisplayName', '$x$ (stable)') %y1
    F_graz
    indF = find(data_all(:,1) < F_graz)
    indF = indF(end)
    
    plot(data_all(ind+1:indF,1),data_all(ind+1:indF,6),'r-', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    plot(data_all(indF+1:end,1),data_all(indF+1:end,6),'r--', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    
    plot(data_all(1:indF,1),data_all(1:indF,3),'k-', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1
    plot(data_all(indF+1:2506,1),data_all(indF+1:2506,3),'k--', 'LineWidth', lw, 'DisplayName', '$y$ (unstable)') %y1

    % plot(data_all(ind+1:end,1),data_all(ind+1:end,6),'r-', 'LineWidth', lw, 'DisplayName', '$x$ (unstable)') %y1
    s = (data_all(ind-1,3)-1)/(data_all(ind-1,3)-data_all(ind,3));
    Fyt_impact_SN = (1-s)*data_all(ind-1,[1,6,7, 13]) + s*data_all(ind,[1,6,7, 13])
    ylim([-.4, 1.4]);
    xlim([0.3628, 0.3631]);
    yline(1, '-', 'LineWidth', lw-2, 'HandleVisibility', 'off', 'Color', [0.5, 0.5, 0.5])
    xlabel('$\mathcal{A}$', 'Interpreter', 'latex', 'FontSize',30);
    % ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 30, 'Rotation', 0);
    title('$\omega=0.799$', 'FontSize', 30, 'Interpreter','latex');
    plot(Fyt_impact_SN(1), 1,'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [144, 103, 167]/256, 'MarkerSize', ms, 'HandleVisibility', 'off')    
    text(0.362815,0.92,'SN', 'FontSize', 20)
    xline(F_graz, 'LineStyle', '-', 'Color', 'black', 'LineWidth', lw, 'DisplayName', 'grazing')
    xline(F_graz, 'LineStyle', '-', 'Color', 'green', 'LineWidth', lw, 'DisplayName', 'grazing')
    
    text(0.36292,0.83,'$\lambda_1$', 'Interpreter', 'latex', 'FontSize',25)
    text(0.363025,-0.267526,'$\lambda_2$', 'Interpreter', 'latex', 'FontSize',25)
    text(0.363053,0.07,'$y_{\rm imp}$', 'Interpreter', 'latex', 'FontSize',25)

elseif runNum==23
    figure(62)
    set(gcf,'Position',[100 100 500 440])
    set(gca,'Position',[0.14,0.14,0.8,0.8])
    hold on
    grid off
    
    % set(groot, 'defaultLegendInterpreter', 'latex');
        
    % plot(data_all(:,1),data_all(:,6),'b-', 'LineWidth', lw, 'DisplayName', '$x$') %y1
    % plot(data_all(:,1),data_all(:,7),'r-', 'LineWidth', lw) %t1
    % plot(data_all(:,1),data_all(:,2),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_1$') %e1
    % plot(data_all(:,1),data_all(:,3),'k-', 'LineWidth', lw, 'DisplayName', '$\lambda_2$') %e2
    % legend('Location', 'southwest');

    %crude approximation to SN from given data
    % ind = find(data_all(:,3) > 1,1);

    indF = find(data_all(:,1) < F_graz);
    indF = indF(end);

    % plot(data_all(:,1),data_all(:,2),'k-', 'LineWidth', lw, 'HandleVisibility', 'off') %e1
    % plot(data_all(:,1),data_all(:,3),'k-', 'LineWidth', lw, 'HandleVisibility', 'off') %e2


    plot(data_all(1:indF,1),data_all(1:indF,2),'k-', 'LineWidth', lw, 'MarkerSize',6) %y1
    plot(data_all(indF+1:end,1),data_all(indF+1:end,2),'k--', 'LineWidth', lw, 'HandleVisibility', 'off') %y1

    plot(data_all(1:indF,1),data_all(1:indF,3),'k-', 'LineWidth', lw, 'MarkerSize',6) %y1
    plot(data_all(indF+1:end,1),data_all(indF+1:end,3),'k--', 'LineWidth', lw, 'HandleVisibility', 'off') %y1


    % plot(data_all(:,1),data_all(:,6),'b-', 'LineWidth', lw, 'HandleVisibility', 'off') %y1
    plot(data_all(1:indF,1),data_all(1:indF,6),'b-', 'LineWidth', lw, 'MarkerSize',6) %y1
    plot(data_all(indF+1:end,1),data_all(indF+1:end,6),'b--', 'LineWidth', lw, 'HandleVisibility', 'off') %y1

    ylim([-.4, 1.4]);
    xlim([0.3628, 0.3631]);
    
end

set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);


% dlmwrite('oneParamData.csv', data_all, 'delimiter', ',', 'precision', '%.16f')
% F_graz
