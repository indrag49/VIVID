function schemPoinMap


hh = 4*4/3;
arrSize = .2/hh;
componentColor = [.4,0,.8];
impactColor = [.1,.5,.7];

xMin = -1.5;
xMax = 1;
yMin = -1;
yMax = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf
hold on

axis([xMin,xMax,yMin,yMax])
set(gcf,'Position',[500,200,840,420])
set(gca,'Position',[0,0,1,1])
axis off


fill([0, xMax, xMax, 0], [yMin, yMin, yMax, yMax], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5)

plot(xMin+[0,.98*(xMax-xMin)],[0,0],'k-','LineWidth',5/hh)
plot([0,0],yMin+[0,.98*(yMax-yMin)],'k-','LineWidth',5/hh)

ARR = getArrowNew([xMin,xMax;0,0]',2,1,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
fill(ARR([2,1,4],1),ARR([2,1,4],2),'k','EdgeColor','none')

ARR = getArrowNew([0,0;yMin,yMax]',2,1,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
fill(ARR([2,1,4],1),ARR([2,1,4],2),'k','EdgeColor','none')

A = [0,1;-1,-.3];
zEq = [-.7;0];
tInt = .83;
p2 = expm(tInt*A)*([.3;0] - zEq) + zEq;

for k = 1:2
    if k == 1
        tt = linspace(0,tInt,100);
    elseif k == 2
        tt = linspace(tInt,2*pi/imag([1,0]*eig(A)),200);
    end
    zz = zeros(2,length(tt));
    for j = 1:length(tt)
        zz(:,j) = expm(tt(j)*A)*([.3;0] - zEq) + zEq;
    end
    if k == 1
        plot(zz(1,:),zz(2,:),'k--','LineWidth',10/hh)
        ARR = getArrowNew(zz',60,10,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
        fill(ARR(:,1),ARR(:,2),'k','EdgeColor','none')
        p3 = zz(:,1);
    elseif k == 2
        plot(zz(1,:),zz(2,:),'k-','LineWidth',10/hh)
        ARR = getArrowNew(zz',50,1,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
        fill(ARR(:,1),ARR(:,2),'k','EdgeColor','none')
        ARR = getArrowNew(zz',170,5,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
        fill(ARR(:,1),ARR(:,2),'k','EdgeColor','none')
        p4 = zz(:,end);
    end
end

A = [0,1.5;-.5,-.1];
zEq = [-.5;0];
for k = 1:2
    if k == 1
        tt = linspace(0,1.28,100);
    elseif k == 2
        tt = linspace(-.9,0,100);
    end
    zz = zeros(2,length(tt));
    for j = 1:length(tt)
        zz(:,j) = expm(tt(j)*A)*([0;-p2(2)] - zEq) + zEq;
    end
    if k == 1
        plot(zz(1,:),zz(2,:),'k--','LineWidth',10/hh)
        ARR = getArrowNew(zz',60,3,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
        fill(ARR(:,1),ARR(:,2),'k','EdgeColor','none')
        p0 = zz(:,end);
        p1 = zz(:,1);
    elseif k == 2
        plot(zz(1,:),zz(2,:),'k-','LineWidth',10/hh)
        ARR = getArrowNew(zz',30,1,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
        fill(ARR(:,1),ARR(:,2),'k','EdgeColor','none')
    end
end

pp = [p0,p1,p2,p3,p4];
plot(pp(1,:),pp(2,:),'ko','MarkerFaceColor','k','MarkerSize',35/hh)


ss = linspace(0,1,100);
zz = p1*(1-ss) + p2*ss + .2*[1;0]*ss.*(1-ss);
plot(zz(1,:),zz(2,:),'-','LineWidth',10/hh,'Color',impactColor)
ARR = getArrowNew(zz',40,1,arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
fill(ARR(:,1),ARR(:,2),impactColor,'EdgeColor','none')

ind = 25;
% pointer line for instantaneous velocity reversal
% plot(zz(1,ind)+[0,.4],zz(2,ind)+[0,.45],'k-','LineWidth',5/hh)
zz = p0*(1-ss) + p3*ss - .25*[0;1]*ss.*(1-ss);
plot(zz(1,:),zz(2,:),'-','LineWidth',5/hh,'Color',componentColor)
ARR = getArrowNew(zz',85,1,.7*arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
fill(ARR(:,1),ARR(:,2),componentColor,'EdgeColor','none')

zz = p3*(1-ss) + p4*ss - .35*[0;1]*ss.*(1-ss);
plot(zz(1,:),zz(2,:),'-','LineWidth',5/hh,'Color',componentColor)
ARR = getArrowNew(zz',80,1,.7*arrSize,xMax-xMin,yMax-yMin,get(gca,'Position'),get(gcf,'Position'),.5,.5,0);
fill(ARR(:,1),ARR(:,2),componentColor,'EdgeColor','none')

text(.86,.055,'$x$','interpreter','latex','FontSize',25)
text(.035,.92,'$y$','interpreter','latex','FontSize',25)

text(.4,-.15,'$P_{{\rm disc},R}$','interpreter','latex','FontSize',25,'Color',componentColor)
text(-.27,-.14,'$P_{\rm global}$','interpreter','latex','FontSize',25,'Color',componentColor)

text(0.01, 0.74, '$(y_1,\mathbf{z}_1)$','interpreter','latex','FontSize',25)
text(0.65, -0.1, '$(x_0,\mathbf{z}_0)$','interpreter','latex','FontSize',25)
text(0.29, 0.1, '$(x_3,\mathbf{z}_3)$','interpreter','latex','FontSize',25)
text(-0.31, 0.1, '$(x_4,\mathbf{z}_4)$','interpreter','latex','FontSize',25)
text(0.01, -0.72, '$(y_2,\mathbf{z}_1)$','interpreter','latex','FontSize',25)

text(0.065, 0.22, '$\Phi$','interpreter','latex','FontSize',25)
text(-0.6, -0.06, '$\Pi$','interpreter','latex','FontSize',25)
text(-0.08, 0.4, '$\Sigma$','interpreter','latex','FontSize',25)

