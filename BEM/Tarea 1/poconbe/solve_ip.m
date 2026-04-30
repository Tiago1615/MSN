function [ipT,ipq]=solve_ip(nx,en,nT,nq,ipx)
% SOLVE_IP Calculate the solution at internal points from the solution 
%          throughout the boundary.
%
%   [ipT,ipq]=SOLVE_IP(nx,en,nT,nq,ipx) Calculate the solution T' and 
%   q' (q'_x and q'_y) at internal points from the nodal solutions T (nT) 
%   and q' (nq) of boundary.
%
%   See also POCONBE, UNMAP, WRITE_RESULTS.

  % Set the number of elements and functional nodes (the same number)
  ne=length(en(1,:));
  % Set the number of internal points
  nip=length(ipx(1,:));
  % Preallocate solution at internal points
  ipT=zeros(nip,1);
  ipq=zeros(2,nip);
  %
  % Loop through elements to integrate
  %
  for kj=1:ne
    %
    % Calculate required element data for integration
    %
    x1=nx(:,en(1,kj));
    x2=nx(:,en(2,kj));
    r12=x2-x1;
    L=sqrt(r12(1)^2+r12(2)^2);
    n(1)=r12(2)/L;
    n(2)=-r12(1)/L;
    %
    % Loop through internal points to collocate
    %
    for ki=1:nip
      % Coordinates of the collocation point
      x_i=ipx(:,ki);
      % Calculate h, g, m and l
      [h,g,m,l]=integrate_exterior_hgml(x_i,x1,x2,L,n);
      % Assemble h, g, m and l
      ipT(ki)=ipT(ki)+g*nq(kj)-h*nT(kj);
      ipq(1,ki)=ipq(1,ki)+l(1)*nq(kj)-m(1)*nT(kj);
      ipq(2,ki)=ipq(2,ki)+l(2)*nq(kj)-m(2)*nT(kj);
    end
  end
  
end
