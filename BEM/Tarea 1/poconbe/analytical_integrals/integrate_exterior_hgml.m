function [h,g,m,l]=integrate_exterior_hgml(x_i,x1,x2,L,n)
% INTEGRATE_EXTERIOR_HGML Calculates kernels for u and q BIE when the collocation point is
%                       not at the element.
%
%   [h,g,m,l]=INTEGRATE_EXTERIOR_HGML(x_i,x1,x2,L,n) Calculates h, g, m and l kernels 
%   for u and q BIE when the collocation point is not at the element. x_i is the 
%   collocation point coordinates, x1 and x2 are the vectors of element nodal 
%   coordinates, L is the element length and n is the element unit normal.
%
%   See also BUILD, SOLVE_IP.

  a=0.5*(x1+x2)-x_i;
  b=0.5*(x2-x1);
  c0=a(1)^2+a(2)^2;
  c1=2*(a(1)*b(1)+a(2)*b(2));
  c2=0.25*L^2;  
  d0=c1^2-4*c0*c2;
  I0=integral_lnr(a,b,c0,c1,c2,d0,L);
  I1=integral_drdxk_div_r(a,b,c0,c1,c2,d0,L);
  I2=integral_1_div_r2(a,b,c0,c1,c2,d0,L);
  I3=integral_drdxj_drdxk_div_r2(a,b,c0,c1,c2,d0,L);
  cte=1/(2*pi);
  g=-cte*I0;
  h=-cte*dot(I1,n);  
  m(1)=cte*(-2*(I3(1,1)*n(1)+I3(1,2)*n(2))+I2*n(1));
  m(2)=cte*(-2*(I3(2,1)*n(1)+I3(2,2)*n(2))+I2*n(2));
  l(1)=cte*I1(1);
  l(2)=cte*I1(2);

end
