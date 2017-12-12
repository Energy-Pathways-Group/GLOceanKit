%--------------------------------------------------
% quick look at flow fields
%--------------------------------------------------
  clear all; close all
  FIGPOS = [405 145 1009 635];
  do_dye=0;

  numprocs = 16;  % number of processors used in computing the solutions
  Ly = 1.6e4;     % [m]
  Lz = 2400;     % [m]
  nu = 2.e-6;    % [m2/s]
  kappa=nu;      % [m2/s]
  dt = 0.015;    % [s]
  rho0 = 1000;   % [kg/m3]
  g = 9.81;      % [m/s2]
 


  g=9.81;rho0=1000;    % [m/s2],[kg/m3]
  dRdz=0.0;            % [kg/m4]

 
 for slice=16400:16400         % which time slice to grab from the netcdf file...

  !rm -f slice_*.nc      % matlab passes this directly to the shell
  cmd=['!./cat_slices.pl ' num2str(slice+1) ' ' num2str(numprocs)];
  eval(cmd);               % pass command to matlab, which passes it to the shell because of ! at beginning

  cnum=sprintf('%.6i',slice);
  fname=['slice_' cnum '.nc'];

  % use snctools in matlab to extract data from netcdf file...
  y=double( nc_varget(fname,'y') );
  z=double( nc_varget(fname,'z') );
  time=double( nc_varget(fname,'time') );
  v=double( nc_varget(fname,'v') );
  w=double( nc_varget(fname,'w') );
  rho=double( nc_varget(fname,'s1') ); 
  if(do_dye==1), s2=double( nc_varget(fname,'s2') ); end

  ny=length(y); dy=y(2)-y(1);     % determine array sizes...
  nz=length(z); dz=z(2)-z(1);     %  ""
  
   

  


  fh=figure(1);
   clf
   set(fh,'Position',FIGPOS)
   %set(fh,'Color',[0 0 0],'InvertHardcopy','off')

   contour(y,z,rho,32)
   colorbar 
   colormap(flipud(jet))
   hold on
    inc=4;
    quiver(y(1:inc:end),z(1:inc:end),v(1:inc:end,1:inc:end),w(1:inc:end,1:inc:end),2,'k')
   hold off

   xlabel('y/\delta','FontName','Times','FontSize',14)
   ylabel('z/\delta','FontName','Times','FontSize',14)
   title('density + velocity arrows','FontName','Times','FontSize',14)
   set(gca,'FontName','Times','FontSize',14)
 

end   % end loop through slices
!rm -f slice_00*.nc        % clean up

figure(1)
print -dpng sample_image.png


