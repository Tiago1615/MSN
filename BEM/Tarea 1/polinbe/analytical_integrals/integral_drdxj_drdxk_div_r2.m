function I=integral_drdxj_drdxk_div_r2(ap,bp,c0p,c1p,d0p,L)
  % Solution to the integral: int_{-1}^{1} r_{,j}·r_{,k}/r^2·L/2·dxi
  if abs(d0p)<1.e-14
    if abs(0.5*c1p)<1+1e-14
      I=[0 0;0 0];
    else    
      ajp=ap(1);
      bjp=bp(1);
      akp=ap(1);
      bkp=bp(1);
      I(1,1)=(64*(4*ajp*akp+12*bjp*bkp-8*(akp*bjp+ajp*bkp)*c1p...
        +(3*ajp*akp+bjp*bkp)*c1p^2))/(3.*(-4+c1p^2)^3*L);
      ajp=ap(1);
      bjp=bp(1);
      akp=ap(2);
      bkp=bp(2);
      I(1,2)=(64*(4*ajp*akp+12*bjp*bkp-8*(akp*bjp+ajp*bkp)*c1p...
        +(3*ajp*akp+bjp*bkp)*c1p^2))/(3.*(-4+c1p^2)^3*L);
      I(2,1)=I(1,2);
      ajp=ap(2);
      bjp=bp(2);
      akp=ap(2);
      bkp=bp(2);
      I(2,2)=(64*(4*ajp*akp+12*bjp*bkp-8*(akp*bjp+ajp*bkp)*c1p...
        +(3*ajp*akp+bjp*bkp)*c1p^2))/(3.*(-4+c1p^2)^3*L);
    end
  else
    dis=sqrt(-d0p);
    den=((1+c0p-c1p)*(1+c0p+c1p)*dis^3*L);
    atandif=atan((-2+c1p)/dis)-atan((2+c1p)/dis);
    ajp=ap(1);
    bjp=bp(1);
    akp=ap(1);
    bkp=bp(1);
    I(1,1)=(-4*(dis*(2*(1+c0p)*(-(ajp*akp)+bjp*bkp*c0p)...
      -(akp*bjp+ajp*bkp)*(-1+c0p)*c1p+(ajp*akp-bjp*bkp)*c1p^2)...
      +(1+c0p-c1p)*(1+c0p+c1p)*(2*ajp*akp+2*bjp*bkp*c0p...
      -(akp*bjp+ajp*bkp)*c1p)*atandif))/den;
    ajp=ap(1);
    bjp=bp(1);
    akp=ap(2);
    bkp=bp(2);
    I(1,2)=(-4*(dis*(2*(1+c0p)*(-(ajp*akp)+bjp*bkp*c0p)...
      -(akp*bjp+ajp*bkp)*(-1+c0p)*c1p+(ajp*akp-bjp*bkp)*c1p^2)...
      +(1+c0p-c1p)*(1+c0p+c1p)*(2*ajp*akp+2*bjp*bkp*c0p...
      -(akp*bjp+ajp*bkp)*c1p)*atandif))/den;
    I(2,1)=I(1,2);
    ajp=ap(2);
    bjp=bp(2);
    akp=ap(2);
    bkp=bp(2);
    I(2,2)=(-4*(dis*(2*(1+c0p)*(-(ajp*akp)+bjp*bkp*c0p)...
      -(akp*bjp+ajp*bkp)*(-1+c0p)*c1p+(ajp*akp-bjp*bkp)*c1p^2)...
      +(1+c0p-c1p)*(1+c0p+c1p)*(2*ajp*akp+2*bjp*bkp*c0p...
      -(akp*bjp+ajp*bkp)*c1p)*atandif))/den;
  end
end

