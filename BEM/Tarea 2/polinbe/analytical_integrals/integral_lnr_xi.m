function I=integral_lnr_xi(ap,bp,c0p,c1p,d0p,L)
  % Solution to the integral: int_{-1}^{1} ln(r)·xi·L/2·dxi 
  if abs(d0p)<1.e-14
    if abs(0.5*c1p-1)<1.e-14
      I= L/2;
    elseif abs(0.5*c1p+1)<1.e-14
      I=-L/2;
    elseif 0.5*c1p>1
      I=(L*(c1p-((-4+c1p^2)*atanh(2/c1p))/2.))/4.;
    elseif 0.5*c1p<-1
      I=(L*(2*c1p-c1p^2*atanh(2/c1p)+log((2+c1p)^2/(-2+c1p)^2)))/8.;
    else
      I=(L*(-(c1p^2*atanh(c1p/2.))+log((-2+c1p)^(-2))+2*(c1p+log(2+c1p))))/8.;
    end
  else
    dis=sqrt(-d0p);
    I=(L*(4*c1p+2*c1p*dis*atan((-2+c1p)/dis)...
      -2*c1p*dis*atan((2+c1p)/dis)...
      +(-2-2*c0p+c1p^2)*log((1+c0p-c1p)/(1+c0p+c1p))))/16.;
  end
end
