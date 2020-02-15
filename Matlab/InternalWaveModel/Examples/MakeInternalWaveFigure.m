%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 15e3;
Ly = 15e3;
Lz = 1300;

Nx = 32;
Ny = 32;
Nz = 33;

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 1; % k=0..Nx/2
l0 = 1; % l=0..Ny/2
j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 0.15; % m/s
sign = 1;

period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);

t = 0;
[u,v,w,rho_prime] = wavemodel.VariableFieldsAtTime(t,'u','v','w','rho_prime');

rho_bar = wavemodel.RhoBarAtDepth(wavemodel.z);
rho3d = wavemodel.RhoBarAtDepth(wavemodel.Z) + rho_prime;

deltaX = wavemodel.x(2)-wavemodel.x(1);
minX = min(wavemodel.x);
maxX = max(wavemodel.x+deltaX);

deltaY = wavemodel.y(2)-wavemodel.y(1);
minY = min(wavemodel.y);
maxY = max(wavemodel.y+deltaY);

minZ = min(wavemodel.z);
maxZ = max(wavemodel.z);

[X,Y,Z] = meshgrid(wavemodel.x,wavemodel.y,wavemodel.z);


% FigureSize = [2 2 24 20];
% fig1 = figure('Units', 'centimeters', 'Position', FigureSize);
% % set(gcf,'PaperPositionMode','auto')
% fig1.PaperUnits = 'points';
% fig1.PaperPosition = FigureSize;
% fig1.PaperSize = [FigureSize(3) FigureSize(4)];

figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% specify the density surface that we want to display, and find their color in the colormap
%
colormap(flipud(jet(128)))
%colormap(jet(128))
coloraxis = linspace(0,1,length(colormap));
deltaRho = (max(rho_bar)-min(rho_bar))/3;
density_surface = [min(rho_bar) min(rho_bar)+deltaRho min(rho_bar)+2*deltaRho];
density_color = interp1(coloraxis,colormap, (density_surface-min(rho_bar))./(max(rho_bar)-min(rho_bar)) );


		
hsurfaces = slice(X,Y,Z,rho3d,[min(wavemodel.x)],[max(wavemodel.y)],[min(wavemodel.z)]);
set(hsurfaces,'FaceColor','interp','EdgeColor','none')
caxis([min(rho_bar) max(rho_bar)])

view(30,10);

for i=1:length(density_surface)
    p = patch(isosurface(X,Y,Z,rho3d,density_surface(i)));
    isonormals(X,Y,Z,rho3d,p)
    alpha(p,0.5)
    set(p,'FaceColor',density_color(i,:),'EdgeColor','none');
end

hold on


xlim([minX maxX])
ylim([minY maxY])
zlim([min(wavemodel.z) max(wavemodel.z)])

lighting gouraud
camlight(30,20)
camlight


hold off

axis off

print('InternalWave.png','-dpng','-r300')