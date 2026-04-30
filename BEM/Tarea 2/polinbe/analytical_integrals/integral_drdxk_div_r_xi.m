function I=integral_drdxk_div_r_xi(ap,bp,c0p,c1p,d0p,L)
  % Solution to the integral: int_{-1}^{1} r_{,k}/r·xi·L/2·dxi
  if abs(d0p)<1.e-14
    if abs(0.5*c1p)<1+1e-14
      I=[0 0];
    else
      I=(-4*ap*c1p+4*bp*(-2+c1p^2))/(-4+c1p^2)+2*(ap-bp*c1p)*atanh(2/c1p);
    end
  else
    dis=sqrt(-d0p);
    I=(4*bp*dis...
      +2*(2*bp*c0p+ap*c1p-bp*c1p^2)*(atan((-2+c1p)/dis)-atan((2+c1p)/dis))...
      +dis*(-2*bp*c1p*atanh(c1p/(1+c0p))+ap*log((1+c0p+c1p)/(1+c0p-c1p))))/(2.*dis);
  end
end
