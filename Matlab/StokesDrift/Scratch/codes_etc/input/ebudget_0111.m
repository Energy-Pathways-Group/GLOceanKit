
% program to compute total potential energy budget using  3d data files from Kraig's code
% derivatives are computed in spectral space.
disp('computation of various terms in the energy budget')

%%%%%%%%%%%%% simulation parameters **********************
%%%%%%%%%% must be edited for each run %%%%%%%%%%%%%%%%%%%

gravity= 9.81;
rho0 = 1027
nu=2.5e-6; % m^2/s  Fickian viscosity
kappa2=nu; % m^2/s  Fickian diffusivity
length_scale=1.e3 % m depth of domain
nu6=(1.5e16)/length_scale^6;   % (corrected) hyperviscosity
kappa6=nu6/10.;
DGRAD=4.18e-4;
dt = 31.4*200/16.   % time elapsed between 3d field saves (seconds)
coriolis=2.e-4;
omega=coriolis/0.85;
wave_period=2.*pi/omega;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' compute energies, dissipation from 3d cfd files')
files0=dir('cfd*_000');  % stores all filenames with cfd* into a structure "files"
files1=dir('cfd*_001');
files2=dir('cfd*_002');
files3=dir('cfd*_003');
files4=dir('cfd*_004');
files5=dir('cfd*_005');
files6=dir('cfd*_006');
files7=dir('cfd*_007');
nfiles=size(files0)
nprocs=8;
%nodata_flag=0;
%if ~isempty(files);

% definition of various energy budget terms

t=zeros(nfiles,1); % time corresponding to saved 3d fields
% definition of big arrays
file0=getfield(files0(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file0); % read it
% define various quantities...
nx=size(x,1);
ny=size(y,1);
locnz=size(zi,1);
nz=nprocs*locnz;
u=zeros(nx,ny,nz);
v=zeros(nx,ny,nz);
w=zeros(nx,ny,nz);
rho=zeros(nx,ny,nz);
z=zeros(nz,1);
z1=zeros(locnz,1);
dy=y(2)-y(1);
dx=x(2)-x(1);
Lx=max(x)+dx;
Ly=max(y)+dy;
%z=linspace(0,Lz-(z1(2)-z1(1)),nz);
%dz=z(2)-z(1);
n1=1;
n2=locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file1=getfield(files1(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file1); % read it
n1=n2+1;
n2=2*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file2=getfield(files2(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file2); % read it
n1=n2+1;
n2=3*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file3=getfield(files3(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file3); % read it
n1=n2+1;
n2=4*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file4=getfield(files4(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file4); % read it
n1=n2+1;
n2=5*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file5=getfield(files5(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file5); % read it
n1=n2+1;
n2=6*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file6=getfield(files6(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file6); % read it
n1=n2+1;
n2=7*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;


file7=getfield(files7(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file7); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;


%define z coordinate
dzi=zi(2)-zi(1);
Lz=1000.
z=linspace(0,Lz-dzi,nz);
dz=z(2)-z(1);


%file=input('filename?=    ','s')   % read in input file

% define various quantities...
r2=zeros(nx,ny,nz);
r2scaled=zeros(nx,ny,nz);


% define wavenumbers for triply periodic domain

dkx=2.*pi/Lx;
dky=2.*pi/Ly;
dkz=2.*pi/Lz;

kx(1)=0.;
for i=2:nx/2
kx(i)=(i-1)*dkx;
kx(nx-i+2) = -kx(i);
end% i
ky(1) = 0.0;
ky(nx/2+1)=dx*nx/2;

for j=2:ny/2
ky(j)=(j-1)*dky;
ky(ny-j+2) = -ky(j);
end% j
ky(1) = 0.0;
ky(ny/2+1) = dky*ny/2;

for k=2:nz/2
kz(k)=(k-1)*dkz;
kz(nz-k+2) = -kz(k);
end% k
kz(1)=0.;
kz(nz/2+1)=dkz*(nz/2);

%%%% maximum wavenumbers (8/9 rule)
kxmax=ceil((sqrt(8./9*(nx/2)^2)+1))*dkx;
kymax=ceil((sqrt(8./9*(ny/2)^2)+1))*dky;
kzmax=ceil((sqrt(8./9*(nz/2)^2)+1))*dkz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:nz
for j=1:ny
for i=1:nx
r2(i,j,k)=kx(i)^2+ky(j)^2+kz(k)^2;
r2scaled(i,j,k)=((kx(i)/kxmax)^2 + (ky(j)/kymax)^2 + (kz(k)/kzmax)^2)^3;
end
end
end

% definition of grid spacings used in hyperdiffusion operator
dxmod=dx*kxmax;
dymod=dy*kymax;
dzmod=dz*kzmax;

%define various arrays
del2rho=zeros(nx,ny,nz); % work array
del6mod=zeros(nx,ny,nz); 
del4mod=zeros(nx,ny,nz);
gradmodx=zeros(nx,ny,nz);
gradmody=zeros(nx,ny,nz);
gradmodz=zeros(nx,ny,nz);
gradx=zeros(nx,ny,nz);
grady=zeros(nx,ny,nz);
gradz=zeros(nx,ny,nz);
% ep energy budget
buoy_flux=zeros(nfiles,1); % buoyancy flux
adv_flux=zeros(nfiles,1);  % advective flux 
ep_surf_adv= zeros(nfiles,1); %ep advective surface flux
ep_surf_dif= zeros(nfiles,1); %ep diffusive surface flux
ep_surf_hdif= zeros(nfiles,1); % ep hyperdiffusive surface flux
phi_i=zeros(nfiles,1);
phih_i=zeros(nfiles,1);
% eb energy budget
diap_flux=zeros(nfiles,1); % diapycnal flux
eb_adv= zeros(nfiles,1); %eb advective flux
eb_dif= zeros(nfiles,1); %eb diffusive flux
eb_hdif= zeros(nfiles,1); %eb hyperdiffusive flux
dissh_spec=zeros(nfiles,1);

ke=zeros(nfiles,1);
eb=zeros(nfiles,1);
ep=zeros(nfiles,1);
epa=zeros(nfiles,1);

rhobar=zeros(nx,ny,nz);  % 3d mean density 
rhoprime=zeros(nx,ny,nz); % perturbation density
zz=zeros(nx,ny,nz); % 3d vertical coordinate

% sorted variables
rhostar=zeros(nx,ny,nz); % sorted density
zstar=zeros(nx,ny,nz,1); % sorted z-coordinate

for k=1:nz
zz(:,:,k)=z(k);
rhobar(:,:,k)=DGRAD*(Lz-zz(:,:,k));
zstar(:,:,k)=z(k);
end;
drhobar=-DGRAD*ones(nx,ny,nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP IN TIME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nfiles;
file0=getfield(files0(n,1),'name') % get the filename associated with first file
n,file0
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file0); % read it
% define various quantities...
n1=1;
n2=locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file1=getfield(files1(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file1); % read it
n1=n2+1;
n2=2*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file2=getfield(files2(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file2); % read it
n1=n2+1;
n2=3*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file3=getfield(files3(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file3); % read it
n1=n2+1;
n2=4*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file4=getfield(files4(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file4); % read it
n1=n2+1;
n2=5*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file5=getfield(files5(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file5); % read it
n1=n2+1;
n2=6*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file6=getfield(files6(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file6); % read it
n1=n2+1;
n2=7*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;


file7=getfield(files7(n,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file7); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

rhoprime(:)=rho(:)-rhobar(:);
t(n)=(n-1)*dt; % in seconds

% compute sorted profiles
rho=reshape(rho,nx*ny*nz,1);
[rhostar,I]=sort(rho,1,'descend'); % sort in descending order
rhostar=reshape(rhostar,nx,ny,nz,1);
zstar=reshape(zstar,nx*ny*nz,1);
zstar(I)=zz; % attributes each index to a particular z position in x,y,z space
zstar=reshape(zstar,nx,ny,nz);
rho=reshape(rho,nx,ny,nz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kinetic energy per unit volume
ke(n) = 0.5*rho0*(sum(u(:).^2+v(:).^2+w(:).^2))/(nx*ny*nz);

% total potential energy per unit volume
ep(n)=gravity*sum((rho(:)).*zz(:))/(nx*ny*nz); 

% background potential energy per unit volume
eb(n) = gravity*sum(rho(:).*zstar(:))/(nx*ny*nz);

% compute the various background potential energy budget terms
% first compute hyperdiffusive term
%temp=6*del2(rho,dxmod,dymod,dzmod); % factor of 6 because it's 3D (3x2)
frhoprime=fftn(rhoprime);
fu=fftn(u);
fv=fftn(v);
fw=fftn(w);
dissh_spec(n) = nu6*sum(-r2scaled(:).*(abs(fu(:)).^2+abs(fv(:)).^2+abs(fw(:)).^2));



for k=1:nz
for j=1:ny
for i=1:nx
del2rho(i,j,k) = -(kx(i)^2+ky(j)^2+kz(k)^2)*frhoprime(i,j,k);
del4mod(i,j,k) =(-(kx(i)/kxmax)^2-(ky(j)/kymax)^2-(kz(k)/kzmax)^2)^2*frhoprime(i,j,k);
del6mod(i,j,k) = (-(kx(i)/kxmax)^2-(ky(j)/kymax)^2-(kz(k)/kzmax)^2)^3*frhoprime(i,j,k); 
end% i
end% j
end% k
%specifying 'symmetric' imposes conjugate symmetry where it may not be present due to round-off error


% compute gradient in normalized coordinate system
for k=1:nz
for j=1:ny
for i=1:nx
gradmodx(i,j,k) = complex(0,1)*(kx(i)/kxmax)*del4mod(i,j,k);
gradmody(i,j,k) = complex(0,1)*(ky(j)/kymax)*del4mod(i,j,k);
gradmodz(i,j,k) = complex(0,1)*(kz(k)/kzmax)*del4mod(i,j,k);
gradx(i,j,k) =  complex(0,1)*kx(i)*frhoprime(i,j,k);
grady(i,j,k) = complex(0,1)*ky(j)*frhoprime(i,j,k);
gradz(i,j,k) = complex(0,1)*kz(k)*frhoprime(i,j,k);
end;
end;
end;
% transform back to physical space 
del4mod=ifftn(del4mod,'symmetric');  
del6mod=ifftn(del6mod,'symmetric');
del2rho=ifftn(del2rho,'symmetric');
gradmodx=ifftn(gradmodx,'symmetric');
gradmody=ifftn(gradmody,'symmetric');
gradmodz=ifftn(gradmodz,'symmetric');
gradx=ifftn(gradx,'symmetric');
grady=ifftn(grady,'symmetric');
gradz=ifftn(gradz,'symmetric');
%specifying 'symmetric' imposes conjugate symmetry where it may not be present due to round-off error


% advective, diffusive and hyperdiffusive surface fluxes
eb_adv(n) = -gravity*sum((zstar(:).*gradx(:).*u(:)))/(nx*ny*nz)...
-gravity*sum(zstar(:).*grady(:).*v(:))/(nx*ny*nz)...
-gravity*sum(zstar(:).*(drhobar(:)+gradz(:)).*w(:))/(nx*ny*nz);

eb_surf_dif(n)= kappa2*gravity*sum(sum(zstar(:,:,nz).*gradz(:,:,nz)-zstar(:,:,1).*gradz(:,:,1)))/(nx*ny*dz*nz)...
+kappa2*gravity*sum(sum(zstar(nx,:,:).*gradz(nx,:,:)-zstar(1,:,:).*gradz(1,:,:)))/(dx*nx*ny*nz)...
+kappa2*gravity*sum(sum(zstar(:,ny,:).*gradz(:,ny,:)-zstar(:,1,:).*gradz(:,1,:)))/(nx*ny*dy*nz);

eb_surf_hdif(n)=kappa6*gravity*(sum(sum(zstar(nx,:,:).*gradmodx(nx,:,:)-zstar(1,:,:).*gradmodx(1,:,:)))/(dxmod*nx*ny*nz))...
+kappa6*gravity*(sum(sum(zstar(:,ny,:).*gradmody(:,ny,:)-zstar(:,1,:).*gradmody(:,1,:)))/(dymod*nx*ny*nz))...
+kappa6*gravity*(sum(sum(zstar(:,:,nz).*gradmodz(:,:,nz)-zstar(:,:,1).*gradmodz(:,:,1)))/(dzmod*nx*ny*nz));

diap_flux(n)= -eb_surf_dif(n)-eb_surf_hdif(n)+kappa6*gravity*sum(zstar(:).*del6mod(:))/(nx*ny*nz)+kappa2*gravity*sum(zstar(:).*del2rho(:))/(nx*ny*nz)

buoy_flux(n) = gravity*sum(rhoprime(:).*w(:))/(nx*ny*nz);

phi_i(n)=-kappa2*gravity*(rhobar(nx/2,ny/2,nz)-rhobar(nx/2,ny/2,1));
phih_i(n)=-kappa6*gravity*(sum(sum(zz(nx,:,:).*gradmodx(nx,:,:)-zz(1,:,:).*gradmodx(1,:,:)))/(dxmod*nx*ny*nz))...
-kappa6*gravity*(sum(sum(zz(:,ny,:).*gradmody(:,ny,:)-zz(:,1,:).*gradmody(:,1,:)))/(dymod*nx*ny*nz))...
-kappa6*gravity*(sum(sum(zz(:,:,nz).*gradmodz(:,:,nz)-zz(:,:,1).*gradmodz(:,:,1)))/(dzmod*nx*ny*nz));
end % n
t_wp=t/wave_period;
dissh_spec = 2*rho0*dissh_spec/(nx*ny*nz)^2; %in spectral space, hence the added factors....also a factor of rho0 needed
clear u v w rho;
save 'ebudget_0111.mat' 


