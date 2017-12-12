
% program to concatenate data files from Kraig's code
disp('concatenating data files')

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
files0=dir('XZplane_023000*-000.nc');  % stores all filenames with cfd* into a structure "files"
files0=dir('XZplane*-000.nc');
files1=dir('XZplane*-001.nc');
files2=dir('XZplane*-002.nc');
files3=dir('XZplane*-003.nc');
files4=dir('XZplane*-004.nc');
files5=dir('XZplane*-005.nc');
files6=dir('XZplane*-006.nc');
files7=dir('XZplane*-007.nc');
files8=dir('XZplane*-008.nc');
files9=dir('XZplane*-009.nc');
files10=dir('XZplane*-010.nc');
files11=dir('XZplane*-011.nc');
files12=dir('XZplane*-012.nc');
files13=dir('XZplane*-013.nc');
files14=dir('XZplane*-014.nc');
files15=dir('XZplane*-015.nc');

nfiles=size(files0)
nprocs=16;
%nodata_flag=0;
%if ~isempty(files);

% definition of various energy budget terms

%t=zeros(nfiles,1); % time corresponding to saved 3d fields
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

file8=getfield(files8(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file8); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file9=getfield(files9(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file9); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;
file9=getfield(files9(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file9); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file10=getfield(files10(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file10); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file11=getfield(files11(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file11); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file12=getfield(files12(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file12); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file13=getfield(files13(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file13); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file14=getfield(files14(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file14); % read it
n1=n2+1;
n2=8*locnz;
u(:,:,n1:n2)=ui;
v(:,:,n1:n2)=vi;
w(:,:,n1:n2)=wi;
rho(:,:,n1:n2)=rhoi;

file15=getfield(files15(1,1),'name') % get the filename associated with first file
[ui,vi,wi,rhoi,x,y,zi]=rnc_kraig(file15); % read it
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


