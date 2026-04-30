function I=integral_1_div_r2(a,b,c0,c1,c2,d0,L)
  % Solution to the integral: int_{-1}^{1} 1/r^2·L/2·dxi
  if (d0>-1.e-15)
    I=4*c2*L/(c1^2-4*c2^2);
  else
    dis=sqrt(-d0);
    I=(L*(-atan((c1-2*c2)/dis)+atan((c1+2*c2)/dis)))/dis;
  end
end
