function s = S(tt, b, omega, F)
s = [b*omega*F*sin(omega*tt)/(omega^4+(b^2-2)*omega^2+1) + F*(1-omega^2)*cos(omega*tt)/(omega^4+(b^2-2)*omega^2+1) - 1; b*omega^2*F*cos(omega*tt)/(omega^4+(b^2-2)*omega^2+1)- F*omega*(1-omega^2)*sin(omega*tt)/(omega^4+(b^2-2)*omega^2+1)];
