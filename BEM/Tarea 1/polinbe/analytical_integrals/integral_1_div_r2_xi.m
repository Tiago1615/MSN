function I=integral_1_div_r2_xi(ap,bp,c0p,c1p,d0p,L)
  % Solution to the integral: int_{-1}^{1} 1/r^2·xi·L/2·dxi
  if abs(d0p)<1.e-14
    if abs(0.5*c1p)<1+1e-14
      I=0;
    else
      I=(4*((-2*c1p)/(-4+c1p^2)+atanh(2/c1p)))/L;
    end
  else
    dis=sqrt(-d0p);
    I=(2*((c1p*(atan((-2+c1p)/dis)-atan((2+c1p)/dis)))/dis...
      +log((1+c0p+c1p)/(1+c0p-c1p))/2.))/L;
  end
end
