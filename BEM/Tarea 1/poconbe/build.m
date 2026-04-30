function [A,b]=build(Ta,mk,tcc,vcc,nx,en,ec,xf,nf)
  % BUILD Build matrix A and the vector b of the linear system of equations.
  %
  %   [A,b]=BUILD(mk,tcc,vcc,nx,en,ecxf,nf) Takes the model parameters and build 
  %   the matrix A and vector b.
  %
  %   See also POCONBE, READ_MODEL, UNMAP.

  % Set the number of elements and functional nodes (the same number)
  ne=length(en(1,:));
  % Preallocate A and b
  A=zeros(ne,ne);
  b=zeros(ne,1);
  %
  % Loop through the elements to integrate
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
    % Loop through the functional nodes to collocate
    %
    for ki=1:ne
      %
      % Calculate h and g
      %
      if ki==kj
        h=1/2;
        g=L/(2*pi)*(1-log(L/2));
      else
        x_i=xf(:,ki);
        [h,g]=integrate_exterior_hg(x_i,x1,x2,L,n);
      end
      %
      % Assemble h and g into A and b according to the boundary condition
      %
      switch tcc(ec(kj))
        %
        % T = T0  <==>  T' = T0 - Ta
        %
        %   T' = T0 - Ta   known
        %   q'             unknown
        %
        case 0
          A(ki,kj)=-g;
          b(ki)=b(ki)-h*(vcc(ec(kj))-Ta);
        %
        % q = - k · dT/dn = Q0  <==>  q' = -Q0 / k
        %
        %   T'             unknown
        %   q' = -Q0 / k   known
        %
        case 1
          A(ki,kj)=h;
          b(ki)=b(ki)+g*(-vcc(ec(kj))/mk);
        %
        % q = hc · (T-Ta)  <==>  q' = - hc/k · T'
        %
        %   T'             unknown
        %   q' = - hc/k · T' known
        %
        case 2
          A(ki,kj)=h-g*(-vcc(ec(kj))/mk);
        %
        % Invalid BC
        %
        otherwise
          error('Invalid boundary condition type')
      end
    end
  end
  
end
