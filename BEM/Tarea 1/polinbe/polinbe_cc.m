function polinbe_cc(infile,outfile)
  % POLINBE_CC 2D BEM solver with linear elements for heat conduction problem.
  %
  %  Governing equation is two-dimensional Laplace's equation: 
  %
  %    d^2 T   d^2 T
  %    ----- + ----- = 0
  %    d x^2   d y^2
  %
  %  Boundary conditions:
  %
  %    0: Prescribed temperature T0:
  %
  %       T = T0
  %
  %    1: Prescribed heat flux Q0:  
  %    
  %                 d T
  %       q = - k · --- = Q0
  %                 d n
  %
  %       where k is the thermal conductivity of the material. It should
  %       be considered that:
  %
  %       Q0 > 0 (heat out)
  %       Q0 < 0 (heat in)
  %
  %    2: Prescribed heat transfer by convection:
  %  
  %       q = hc · (T-Ta)
  %
  %       where hc is the convective heat transfer coefficient, and Ta is
  %       the ambient temperature.
  %
  %  Units: coherent system for all magnitudes, e.g. SI. 
  %
  %   POLINBE_cc(infile,outfile) reads the model from file "infile" and write
  %   the results in file "outfile".

  % Library of boundary element integrals
  addpath('analytical_integrals');

  % Read model from input file
  [Ta,mk,tcc,vcc,nx,en,ec,ipx]=read_model(infile);
  
  % Plot mesh
  plot_mesh(tcc,nx,en,ec,ipx)
  
  % The following change of variables is used for solving:
  % T' = T-Ta
  % q' = -q/k
  
  % Build the linear system of equations
  [A,b]=build(Ta,mk,tcc,vcc,nx,en,ec);
  
  % Solve the linear system of equations
  x=linsolve(A,b);
  
  % Recover nodal solutions T' and q' from x
  [nT,nq]=unmap(x,Ta,mk,tcc,vcc,nx,en,ec);
  
  % Calculate solution at internal points
  if ~isempty(ipx)
    [ipT,ipq]=solve_ip(nx,en,nT,nq,ipx);
  else
    ipT=[];
    ipq=[];
  end
  
  % Write functional node and internal point results
  write_results(outfile,Ta,mk,nx,nT,nq,ipx,ipT,ipq);  
  
end
