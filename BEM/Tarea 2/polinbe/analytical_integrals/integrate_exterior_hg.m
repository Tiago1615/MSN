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
  
  % Parameters
  a=0.5*(x1+x2)-x_i;
  b=0.5*(x2-x1);
  c0=a(1)^2+a(2)^2;
  c1=2*(a(1)*b(1)+a(2)*b(2));
  c2=0.25*L^2;  
  d0=c1^2-4*c0*c2;
  % Dimensionless parameters
  ap=a/sqrt(c2);
  bp=b/sqrt(c2);
  c0p=c0/c2;
  c1p=c1/c2;  
  d0p=d0/c2^2;
  % Integrals
  I1=integral_lnr(ap,bp,c0p,c1p,d0p,L);
  I2=integral_lnr_xi(ap,bp,c0p,c1p,d0p,L);
  I3=integral_drdxk_div_r(ap,bp,c0p,c1p,d0p,L);
  I4=integral_drdxk_div_r_xi(ap,bp,c0p,c1p,d0p,L);
  cte=1/(2*pi);
  g=-cte*0.5*[I1-I2 I1+I2];
  h=-cte*0.5*[dot(I3,n)-dot(I4,n) dot(I3,n)+dot(I4,n)];
  
end
