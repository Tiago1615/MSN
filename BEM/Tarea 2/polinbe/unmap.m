function [nT,nq]=unmap(x,Ta,mk,tcc,vcc,nx,en,ec)
  % UNMAP Unmap the solution from the linear system of equations and build
  %       the vectors of nodal u and q solutions.
  %
  %   [nT,nq]=UNMAP(x,mk,tcc,vcc,nx,en,ec) Unmap the solution x from the linear 
  %   system of equations and builds the vectors of nodal T' (nT) and q' (nq) 
  %   solutions.
  %
  %   See also POCONBE, BUILD, SOLVE_IP.

  % Set the number of elements and nodes
  ne=length(en(1,:));
  nn=length(nx(1,:));
  % Preallocate T' and q' at functional nodes
  nT=zeros(nn,1);
  nq=zeros(nn,1);
  % Node flag
  node_flag=false(1,nn);
  %
  % Loop through nodes and recover T' and q' from x according to the BC type
  %
  for ke=1:ne
    for kn=1:2
      kj=en(kn,ke);
      if ~node_flag(kj)
        node_flag(kj)=true;
        switch tcc(ec(ke))
          %
          % T' known, q' unknown
          %
          case 0
            nT(kj)=vcc(ec(ke))-Ta;
            nq(kj)=x(kj);
          %
          % T' unknown, q' known
          %
          case 1
            nT(kj)=x(kj);
            nq(kj)=-vcc(ec(ke))/mk;
          %
          % T' unknown, q' = - hc/k · T' known
          %
          case 2
            nT(kj)=x(kj);
            nq(kj)=-vcc(ec(ke))/mk*nT(kj);
          %
          % Invalid BC
          %
          otherwise
            error('Invalid boundary condition type')
        end
      end
    end
  end
  
end
