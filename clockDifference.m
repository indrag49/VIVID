% input a and b taking values modulo T
% return the value of b-a modulo T that is closest to 0
function d = clockDifference(a,b,T)
a = mod(a,T);
b = mod(b,T);
d = b-a;
if d < -T/2
    d = d + T;
elseif d > T/2
    d = d - T;
end
