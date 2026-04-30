function [A,b]=build(Ta,mk,tcc,vcc,nx,en,ec)
  % BUILD Build matrix A and the vector b of the linear system of equations.
  %
  %   [A,b]=BUILD(mk,tcc,vcc,nx,en,ecxf,nf) Takes the model parameters and build 
  %   the matrix A and vector b.
  %
  %   See also POCONBE, READ_MODEL, UNMAP.

  % Set the number of elements and nodes
  ne=length(en(1,:));
  nn=length(nx(1,:));
  % Preallocate A and b
  A=zeros(nn,nn);
  b=zeros(nn,1);
  %
  % Loop through the elements to integrate
  %
  for kj=1:ne
    %
    % Calculate required element data for integration
    %
    kj1=en(1,kj);
    kj2=en(2,kj);
    x1=nx(:,kj1);
    x2=nx(:,kj2);
    r12=x2-x1;
    L=sqrt(r12(1)^2+r12(2)^2);
    n(1)=r12(2)/L;
    n(2)=-r12(1)/L;
    %
    % Loop through the nodes to collocate
    %
    for ki=1:nn
      x_i=nx(:,ki);
      [h,g]=integrate_exterior_hg(x_i,x1,x2,L,n);
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
          A(ki,kj1)=A(ki,kj1)-g(1);
          A(ki,kj2)=A(ki,kj2)-g(2);
          b(ki)=b(ki)-h(1)*(vcc(ec(kj))-Ta);
          b(ki)=b(ki)-h(2)*(vcc(ec(kj))-Ta);
        %
        % q = - k · dT/dn = Q0  <==>  q' = -Q0 / k
        %
        %   T'             unknown
        %   q' = -Q0 / k   known
        %
        case 1
          A(ki,kj1)=A(ki,kj1)+h(1);
          A(ki,kj2)=A(ki,kj2)+h(2);
          b(ki)=b(ki)+g(1)*(-vcc(ec(kj))/mk);
          b(ki)=b(ki)+g(2)*(-vcc(ec(kj))/mk);
        %
        % q = hc · (T-Ta)  <==>  q' = - hc/k · T'
        %
        %   T'             unknown
        %   q' = - hc/k · T' known
        %
        case 2
          A(ki,kj1)=A(ki,kj1)+h(1)-g(1)*(-vcc(ec(kj))/mk);
          A(ki,kj2)=A(ki,kj2)+h(2)-g(2)*(-vcc(ec(kj))/mk);
        %
        % Invalid BC
        %
        otherwise
          error('Invalid boundary condition type')
      end
    end
  end
  %
  % Add free-terms
  %
  %
  % Loop through the nodes
  %
  for ki=1:nn
    % Find the elements connected to the node
    [~,eln] = find(en==ki);
    neln=length(eln);
    % Two elements connected (usual situation)
    if neln==2
      if en(2,eln(1)) == ki
        if en(1,eln(2)) == ki
          r1=nx(:,en(2,eln(1)))-nx(:,en(1,eln(1)));
          r2=nx(:,en(2,eln(2)))-nx(:,en(1,eln(2)));
        else
          error('Some element connected to node %d have incorrect orientation',ki)
        end
      else
        if en(2,eln(2)) == ki
          r1=nx(:,en(2,eln(2)))-nx(:,en(1,eln(2)));
          r2=nx(:,en(2,eln(1)))-nx(:,en(1,eln(1)));
        else
          error('Some element connected to node %d have incorrect orientation',ki)
        end
      end
    % One element connected
    elseif neln==1
      x_i=nx(:,ki);
      % Find another node at the same position
      ki2=0;
      for kj=1:nn
        if kj~=ki && norm(nx(:,kj)-x_i)<1e-12
          ki2=kj;
          break
        end
      end
      if ki2==0
        error('Node %d must share position with another node',ki)
      end
      [~,eln2] = find(en==ki2);
      neln2=length(eln2);
      if neln2==1
        if en(2,eln) == ki
          if en(1,eln2) == ki2
            r1=nx(:,en(2,eln))-nx(:,en(1,eln));
            r2=nx(:,en(2,eln2))-nx(:,en(1,eln2));
          else
            error('Elements connected to nodes %d and %d have incorrect orientation',ki,ki2)
          end
        else
          if en(2,eln2) == ki2
            r1=nx(:,en(2,eln2))-nx(:,en(1,eln2));
            r2=nx(:,en(2,eln))-nx(:,en(1,eln));
          else
            error('Elements connected to nodes %d and %d have incorrect orientation',ki,ki2)
          end
        end
      else
        error('Nodes %d and %d must be connected to only element',ki,ki2)
      end
    else
      error('Node %d is not connected to any element',ki)
    end
    % Calculate free-term from vectors at the vertex
    tc1=-r1/norm(r1);
    tc2=r2/norm(r2);
    n1=[r1(2) -r1(1)]/norm(r1);
    n2=[r2(2) -r2(1)]/norm(r2);
    ns=n1+n2;
    if norm(ns)<1e-12
      error('Elements connected to node %d have incorrect orientation',ki)
    end
    ns=ns/norm(ns);
    alpha=acos(dot(tc1,ns));
    beta=acos(dot(tc2,ns));
    h=1-(alpha+beta)/(2*pi);
    %
    % Assemble h and g into A and b according to the boundary condition
    %
    switch tcc(ec(eln(1)))
      %
      % T = T0  <==>  T' = T0 - Ta
      %
      %   T' = T0 - Ta   known
      %   q'             unknown
      %
      case 0
        b(ki)=b(ki)-h*(vcc(ec(eln(1)))-Ta);
      %
      % q = - k · dT/dn = Q0  <==>  q' = -Q0 / k
      %
      %   T'             unknown
      %   q' = -Q0 / k   known
      %
      case 1
        A(ki,ki)=A(ki,ki)+h;
      %
      % q = hc · (T-Ta)  <==>  q' = - hc/k · T'
      %
      %   T'             unknown
      %   q' = - hc/k · T' known
      %
      case 2
        A(ki,ki)=A(ki,ki)+h;
      %
      % Invalid BC
      %
      otherwise
        error('Invalid boundary condition type')
    end
  end  
  
end
