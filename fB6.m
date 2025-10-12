function [val,NewtonData,flg] = fB6(y,parVec3)

NewtonData = getFP6(y,parVec3);    % [F,e1,e2,x0,t0,y1,t1,y2,t2,x3,t3]
if isempty(NewtonData)
    val = 0;
    NewtonData = [];
    flg = 0;
else
    flg = 1;
    e = NewtonData(2:3);
    critEig = parVec3(7);
    d = real(e) - critEig;
    val = d(abs(d) == min(abs(d)));
    val = val(1);
end
