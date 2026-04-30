function I=integral_lnr(a,b,c0,c1,c2,d0,L)
  % Solution to the integral: int_{-1}^{1} ln(r)·L/2·dxi  
  if (d0>-1.e-15)
    I=L/2*(-2+(c1*atanh((2*c2)/c1))/c2+log(c1^2/(4*c2)-c2));
  else
    dis=sqrt(-d0);
    I=L/(8*c2)*(-8*c2+2*dis*(-atan((c1-2*c2)/dis)+atan((c1+2*c2)/dis))-(c1-2*c2)*log(c0-c1+c2)+(c1+2*c2)*log(c0+c1+c2));
  end
end
