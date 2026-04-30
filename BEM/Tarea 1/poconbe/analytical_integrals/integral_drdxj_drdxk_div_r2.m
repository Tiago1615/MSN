function I=integral_drdxj_drdxk_div_r2(a,b,c0,c1,c2,d0,L)
  % Solution to the integral: int_{-1}^{1} r_{,j}·r_{,k}/r^2·L/2·dxi
  if (d0>-1.e-15)
    aj=a(1);
    bj=b(1);
    ak=a(1);
    bk=b(1);
    I(1,1)=(16*c2^2*((3*aj*ak+bj*bk)*c1^2-8*(ak*bj+aj*bk)*c1*c2+4*(aj*ak+3*bj*bk)*c2^2)*L)/(3*(c1^2-4*c2^2)^3);
    aj=a(1);
    bj=b(1);
    ak=a(2);
    bk=b(2);
    I(1,2)=(16*c2^2*((3*aj*ak+bj*bk)*c1^2-8*(ak*bj+aj*bk)*c1*c2+4*(aj*ak+3*bj*bk)*c2^2)*L)/(3*(c1^2-4*c2^2)^3);
    I(2,1)=I(1,2);
    aj=a(2);
    bj=b(2);
    ak=a(2);
    bk=b(2);
    I(2,2)=(16*c2^2*((3*aj*ak+bj*bk)*c1^2-8*(ak*bj+aj*bk)*c1*c2+4*(aj*ak+3*bj*bk)*c2^2)*L)/(3*(c1^2-4*c2^2)^3);
  else
    dis=sqrt(-d0);
    atandif=(atan((c1-2*c2)/dis)-atan((c1+2*c2)/dis));
    den=(c0-c1+c2)*(c0+c1+c2)*dis^3;
    aj=a(1);
    bj=b(1);
    ak=a(1);
    bk=b(1);
    I(1,1)=(L*(dis*(bj*(ak*c1*(c0-c2)+bk*(c1^2-2*c0*(c0+c2)))+aj*(bk*c1*(c0-c2)+ak*(-c1^2+2*c2*(c0+c2))))...
               -(c0-c1+c2)*(c0+c1+c2)*(2*bj*bk*c0-ak*bj*c1-aj*bk*c1+2*aj*ak*c2)*atandif))/den;
    aj=a(1);
    bj=b(1);
    ak=a(2);
    bk=b(2);
    I(1,2)=(L*(dis*(bj*(ak*c1*(c0-c2)+bk*(c1^2-2*c0*(c0+c2)))+aj*(bk*c1*(c0-c2)+ak*(-c1^2+2*c2*(c0+c2))))...
               -(c0-c1+c2)*(c0+c1+c2)*(2*bj*bk*c0-ak*bj*c1-aj*bk*c1+2*aj*ak*c2)*atandif))/den;
    I(2,1)=I(1,2);
    aj=a(2);
    bj=b(2);
    ak=a(2);
    bk=b(2);
    I(2,2)=(L*(dis*(bj*(ak*c1*(c0-c2)+bk*(c1^2-2*c0*(c0+c2)))+aj*(bk*c1*(c0-c2)+ak*(-c1^2+2*c2*(c0+c2))))...
               -(c0-c1+c2)*(c0+c1+c2)*(2*bj*bk*c0-ak*bj*c1-aj*bk*c1+2*aj*ak*c2)*atandif))/den;
  end
end

