function write_results(outfile,Ta,mk,nx,nT,nq,ipx,ipT,ipq)
% WRITE_RESULTS Write the analysis results to a file.
%
%   WRITE_RESULTS(outfile,xf,nT,nq,ipx,ipT,ipq) Write the analysis results to 
%   the file "outfile".
%
%   See also POCONBE, BUILD, SOLVE_IP.

  % Open file
  f=fopen(outfile,'w');
  % Undo the change of variables
  nT=nT+Ta;
  ipT=ipT+Ta;
  nq=-mk*nq;
  ipq=-mk*ipq;
  % Write nodal solutions
  nn=length(nx(1,:));
  fprintf(f,'%%_node_ ______x1_____ ______x2_____ ______T______ ______q______\n');
  for kj=1:nn
    fprintf(f,' %6d %+.6e %+.6e %+.6e %+.6e\n',kj,nx(1,kj),nx(2,kj),nT(kj),nq(kj));
  end
  fprintf(f,'\n');
  if ~isempty(ipx)
    % Write internal point solutions
    nip=length(ipx(1,:));
    fprintf(f,'%%_intp_ ______x1_____ ______x2_____ ______T______ _____qx1_____ _____qx2_____\n');
    for kj=1:nip
      fprintf(f,' %6d %+.6e %+.6e %+.6e %+.6e %+.6e\n',kj,ipx(1,kj),ipx(2,kj),ipT(kj),ipq(1,kj),ipq(2,kj));
    end
    fprintf(f,'\n');
  end
  % Close file
  fclose(f);
  
end
