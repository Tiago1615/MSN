function I=integral_drdxk_div_r(a,b,c0,c1,c2,d0,L)
  % Solution to the integral: int_{-1}^{1} r_{,k}/r·L/2·dxi
  if (d0>-1.e-15)    
    I=L*(2*(-(b*c1)+2*a*c2)/(c1^2-4*c2^2)+b/c2*atanh(2*c2/c1));
  else
    dis=sqrt(-d0);
    I=L*(2*(b*c1-2*a*c2)*(atan((c1-2*c2)/dis)-atan((c1+2*c2)/dis))+b*dis*(-log(c0-c1+c2)+log(c0+c1+c2)))/(4*c2*dis);
  end
end
