function [Ta,mk,tcc,vcc,nx,en,ec,ipx,xf,nf]=read_model(infile)
  % READ_MODEL Read the model input file and return the model variables
  %
  %   [mk,tcc,vcc,nx,en,ec,ipx]=READ_MODEL(infile) Reads the model from file 
  %   "infile" and returns the ambient temperature (Ta), thermal conductivity 
  %   (mk), vector of boundary condition types (tcc) and values (vcc) by 
  %   boundary, the coordinates of each node (nx), the nodes of each element %
  %   (en), the element's boundary (ec), the coordinates of each internal point 
  %   (xf), and the coordinate of each boundary functional node (xf).
  %
  %   The file format is the following:
  %     <ambient temperature>
  %     <material thermal conductivity>
  %     <number of boundaries>
  %     <Boundary 1 condition: 0, 1 or 2> <BC value: T0, Q0 or hc>
  %     (...) as many lines as the number of boundaries
  %     <number of nodes>
  %     <x1 coordinate of node 1> <x2 coordinate of node 1>
  %     (...) as many lines as the number of nodes
  %     <number of elements>
  %     <element 1 node 1> <element 1 node 2> <element 1 boundary>
  %     (...) as many lines as the number of nodes
  %     <number of internal points>
  %     <x1 coordinate of internal point 1> <x2 coordinate of internal point 1>
  %     (...) as many lines as the number of nodes
  %
  %   See also POCONBE, BUILD

  %
  % Read file
  %
  % Open file
  f=fopen(infile,'r');
  % Read the ambient temperature
  Ta=fscanf(f,'%g',1);
  % Read the thermal conductivity of the material
  mk=fscanf(f,'%g',1);
  % Read the number of boundaries
  nb=fscanf(f,'%d',1);
  % Read boundary condition type and value
  M=fscanf(f,'%f',[2,nb]);
  tcc=int32(M(1,:));
  vcc=M(2,:);
  % Read the number of nodes
  nn=fscanf(f,'%d',1);
  % Read the node coordinates
  nx=fscanf(f,'%f',[2,nn]);
  % Read the number of elements
  ne=fscanf(f,'%d',1);
  % Read the element nodes and boundary which it belongs
  M=fscanf(f,'%d',[3,ne]);
  en=M(1:2,:);
  ec=M(3,:);
  % Read the number of internal points
  nip=fscanf(f,'%d',1);
  if nip>0
    % Read the internal point coordinates 
    ipx=fscanf(f,'%f',[2,nip]);
  else
    ipx=[];
  end
  % Close file
  fclose(f);
  %
  % Some additional calculations
  %
  % Coordinates and unit normal of each functional node at each element
  xf=zeros(2,ne);
  nf=zeros(2,ne);
  for kj=1:ne
    x1=nx(:,en(1,kj));
    x2=nx(:,en(2,kj));
    xf(:,kj)=(x1+x2)/2;
    r12=x2-x1;
    L=sqrt(r12(1)^2+r12(2)^2);
    nf(1,kj)=r12(2)/L;
    nf(2,kj)=-r12(1)/L;
  end

end
