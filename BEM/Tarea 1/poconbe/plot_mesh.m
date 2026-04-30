function plot_mesh(tcc,nx,en,ec,ipx,xf)
%PLOT_MESH Plot mesh
%   Plot of nodes and element by color for each boundary

  figure
  axis equal
  hold on
  CM=jet(length(tcc));
  for kn=1:length(nx(1,:))
    if kn==1
      hp1=plot(nx(1,kn),nx(2,kn),'ko','DisplayName','Geometrical node');
    else
      plot(nx(1,kn),nx(2,kn),'ko','MarkerFaceColor','k')
    end
  end
  for kc=1:length(tcc)
    ke_set=find(ec==kc);
    for ke=ke_set
      x=nx(1,en(1,ke));
      y=nx(2,en(1,ke));
      u=nx(1,en(2,ke))-nx(1,en(1,ke));
      v=nx(2,en(2,ke))-nx(2,en(1,ke));
      if ke==ke_set(1)
        hp2(kc)=quiver(x,y,u,v,0,'MaxHeadSize',0.4,...
          'DisplayName',sprintf('Element (boundary %d)',kc),...
          'color',CM(kc,:));
        hp3(kc)=plot(xf(1,ke),xf(2,ke),'x',...
          'MarkerEdgeColor',CM(kc,:),...
          'DisplayName',sprintf('Functional node (boundary %d)',kc));
      else
        quiver(x,y,u,v,0,'MaxHeadSize',0.4,'color',CM(kc,:))
        plot(xf(1,ke),xf(2,ke),'x','MarkerEdgeColor',CM(kc,:))
      end    
    end
  end
  if ~isempty(ipx)
    for kip=1:length(ipx(1,:))
      if kip==1
        hp4=plot(ipx(1,kip),ipx(2,kip),'ro',...
          'MarkerFaceColor','r','DisplayName','Internal point');
      else
        plot(ipx(1,kip),ipx(2,kip),'ro','MarkerFaceColor','r')
      end 
    end
    legend([hp1 hp2 hp3 hp4])
  else
    legend([hp1 hp2 hp3])
  end
end

