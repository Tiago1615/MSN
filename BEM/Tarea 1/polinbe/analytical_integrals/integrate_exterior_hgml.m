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
  I5=integral_1_div_r2(ap,bp,c0p,c1p,d0p,L);
  I6=integral_1_div_r2_xi(ap,bp,c0p,c1p,d0p,L);
  I7=integral_drdxj_drdxk_div_r2(ap,bp,c0p,c1p,d0p,L);
  I8=integral_drdxj_drdxk_div_r2_xi(ap,bp,c0p,c1p,d0p,L);
  cte=1/(2*pi);
  g=-cte*0.5*[I1-I2 I1+I2];
  h=-cte*0.5*[dot(I3,n)-dot(I4,n) dot(I3,n)+dot(I4,n)];
  % m(j,node)
  m(1,1)=-2*((I7(1,1)-I8(1,1))*n(1)+(I7(1,2)-I8(1,2))*n(2))+(I5-I6)*n(1);
  m(2,1)=-2*((I7(2,1)-I8(2,1))*n(1)+(I7(2,2)-I8(2,2))*n(2))+(I5-I6)*n(2);
  m(1,2)=-2*((I7(1,1)+I8(1,1))*n(1)+(I7(1,2)+I8(1,2))*n(2))+(I5+I6)*n(1);
  m(2,2)=-2*((I7(2,1)+I8(2,1))*n(1)+(I7(2,2)+I8(2,2))*n(2))+(I5+I6)*n(2);
  m=cte*0.5*m;
  % l(j,node) 
  l(1,1)=I3(1)-I4(1);
  l(2,1)=I3(2)-I4(2);
  l(1,2)=I3(1)+I4(1);
  l(2,2)=I3(2)+I4(2);
  l=cte*0.5*l;
  
end
