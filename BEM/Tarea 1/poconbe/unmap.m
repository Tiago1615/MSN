function [nT,nq]=unmap(x,Ta,mk,tcc,vcc,nx,en,ec)
  % UNMAP Unmap the solution from the linear system of equations and build
  %       the vectors of nodal u and q solutions.
  %
  %   [nT,nq]=UNMAP(x,mk,tcc,vcc,nx,en,ec) Unmap the solution x from the linear 
  %   system of equations and builds the vectors of nodal T' (nT) and q' (nq) 
  %   solutions.
  %
  %   See also POCONBE, BUILD, SOLVE_IP.

  % Set the number of elements and functional nodes (the same number)
  ne=length(en(1,:));
  % Preallocate T' and q' at functional nodes
  nT=zeros(ne,1);
  nq=zeros(ne,1);
  %
  % Loop through nodes and recover T' and q' from x according to the BC type
  %
  for kj=1:ne
    switch tcc(ec(kj))
      %
      % T' known, q' unknown
      %
      case 0
        nT(kj)=vcc(ec(kj))-Ta;
        nq(kj)=x(kj);
      %
      % T' unknown, q' known
      %
      case 1
        nT(kj)=x(kj);
        nq(kj)=-vcc(ec(kj))/mk;
      %
      % T' unknown, q' = - hc/k · T' known
      %
      case 2
        nT(kj)=x(kj);
        nq(kj)=-vcc(ec(kj))/mk*nT(kj);
      %
      % Invalid BC
      %
      otherwise
        error('Invalid boundary condition type')
    end
  end
  
end
