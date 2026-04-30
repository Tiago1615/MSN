function I=integral_lnr(ap,bp,c0p,c1p,d0p,L)
  % Solution to the integral: int_{-1}^{1} ln(r)·L/2·dxi
  if abs(d0p)<1.e-14
    if abs(0.5*c1p-1)<1e-14 || abs(0.5*c1p+1)<1e-14
      I=L*(-1+log(L));
    else
      I=(L*(-4+c1p*atanh((4*c1p)/(4+c1p^2))+log((-4+c1p^2)^2/256.)+4*log(L)))/4.;
    end
  else
    dis=sqrt(-d0p);
    I=(L*(-4+dis*(atan((2-c1p)/dis)+atan((2+c1p)/dis))...
      +c1p*atanh(c1p/(1+c0p))+log(((1+c0p-c1p)*(1+c0p+c1p))/16.)+4*log(L)))/4.;
  end
end
