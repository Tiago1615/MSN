function [h,g]=integrate_exterior_hg(x_i,x1,x2,L,n)
% INTEGRATE_EXTERIOR_HG Calculates kernels for u BIE when the collocation point is
%                       not at the element.
%
%   [h,g]=INTEGRATE_EXTERIOR_HG(x_i,x1,x2,L,n) Calculates h and g kernels 
%   for u BIE when the collocation point is not at the element. x_i is the 
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
  I1=integral_lnr(a,b,c0,c1,c2,d0,L);
  I2=integral_drdxk_div_r(a,b,c0,c1,c2,d0,L);
  cte=1/(2*pi);
  g=-cte*I1;
  h=-cte*dot(I2,n);
  
end
