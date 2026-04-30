function I=integral_drdxk_div_r(ap,bp,c0p,c1p,d0p,L)
  % Solution to the integral: int_{-1}^{1} r_{,k}/r·L/2·dxi
  if abs(d0p)<1.e-14
    if abs(0.5*c1p)<1+1e-14
      I=[0 0];
    else
      I=(8*ap-4*bp*c1p)/(-4+c1p^2)+2*bp*atanh(2/c1p);
    end
  else
    dis=sqrt(-d0p);
    I=((-4*ap+2*bp*c1p)*(atan((-2+c1p)/dis)-atan((2+c1p)/dis))...
      +bp*dis*log((1+c0p+c1p)/(1+c0p-c1p)))/(2.*dis);
  end
end
