Lx = 30e3;
Ly = 15e3;
Lz = 5000;

Nx = 256;
Ny = 128;
Nz = 257;

% Nx = 64;
% Ny = 64;
% Nz = 65;
% 
% Nx = 128;
% Ny = 64;
% Nz = 129;

latitude = 31;
N0 = 5.2e-3/2;
t = 12*3600; 1*86400;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

Lh = Lx/32;
Lv = Lz/8;
x0 = Lx/2;
y0 = Ly/2;
z0 = -Lz/2;

X = wavemodel.X;
Y = wavemodel.Y;
Z = wavemodel.Z;

% zeta0 = 100*exp( -((X-x0).^2 + (Y-y0).^2)/(Lh)^2  - ((Z-z0).^2)/(Lv)^2 );
% zeta0 = 100*exp( -((X-x0).^2 + (Y-y0).^2)/(Lh)^2  - ((Z-z0).^2)/(Lv)^2 ).*sin(X/(Lh/4));
% zeta0 = 100*exp( -((X-x0).^2 + (Y-y0).^2)/(Lh)^2  - ((Z-z0).^2)/(Lv)^2 ).*sin(X/(Lh/4)+Z/(Lv/4));
zeta0 = 100*exp( -((X-x0).^2 + (Y-y0).^2)/(Lh)^2  - ((Z-z0).^2)/(Lv)^2 ).*sin(X/(Lh/8)+Z/(Lv/4));


wavemodel.InitializeWithIsopycnalDisplacementField(zeta0);

[u, v, w, rho_prime, zeta, p_wave]= wavemodel.VariableFieldsAtTime(13*3600, 'u', 'v', 'w', 'rho_prime', 'zeta', 'p');

maxU = max(max(max(abs(u))));
maxV = max(max(max(abs(v))));
maxW = max(max(max(abs(w))));
fprintf('Maximum fluid velocity (u,v,w)=(%.2f,%.2f,%.2f) cm/s\n',100*maxU,100*maxV,100*maxW);

dispvar = zeta;
figure
subplot(2,1,1)
pcolor(wavemodel.x,wavemodel.z,squeeze(dispvar(:,Ny/2,:))'),shading flat
subplot(2,1,2)
pcolor(wavemodel.x,wavemodel.y,squeeze(dispvar(:,:,floor(Nz/2)))'),shading flat, axis equal
% figure
% plot(wavemodel.x,squeeze(dispvar(:,Ny/2,floor(Nz/2))));
% 
% figure
% p1 = patch(isosurface(wavemodel.x/1000,wavemodel.y/1000,wavemodel.z,p_wave,0.2));
% isonormals(wavemodel.x/1000,wavemodel.y/1000,wavemodel.z,p_wave,p1)
% alpha(p1,0.5)
% set(p1,'FaceColor',0.8*[1 1 1],'EdgeColor','none');

% figure
% p_mag = 0.1;
% 
% hold on
% [fo,vo] = isosurface(wavemodel.x/1000,wavemodel.y/1000,wavemodel.z,p_wave,p_mag);
% p1 = patch('Faces', fo, 'Vertices', vo);
% p1.FaceColor = 'red';
% p1.EdgeColor = 'none';
% camlight(40,40)                                % create two lights 
% camlight(-20,-10)
% lighting gouraud
% hold on
% [fo,vo] = isosurface(wavemodel.x/1000,wavemodel.y/1000,wavemodel.z,p_wave,-p_mag);
% p2 = patch('Faces', fo, 'Vertices', vo);
% p2.FaceColor = 'blue';
% p2.EdgeColor = 'none';