function ARR = getArrowNew(DATA,j,d,w,xRange,yRange,axPosition,figPosition,ca1,ca2,showConstruct)

YFigRange = 1;
XFigRange = figPosition(3)/figPosition(4)*YFigRange;
XAxRange = axPosition(3)*XFigRange;
YAxRange = axPosition(4)*YFigRange;

A = diag([XAxRange/xRange,YAxRange/yRange]);

vb = (DATA(j-d,:)-DATA(j,:))';
XYb = A*vb;
vb = w*vb/norm(XYb);

XYc = [-XYb(2);XYb(1)];
vc = A\XYc;
vc = w*vc/norm(XYc);

p = DATA(j,:)';
ARR = [p,p+(1+ca1)*vb+ca2*vc,p+vb,p+(1+ca1)*vb-ca2*vc]';

if showConstruct == 1
    plot(p(1)+[0,vb(1);0,vc(1)]',p(2)+[0,vb(2);0,vc(2)]','c*-')
end
