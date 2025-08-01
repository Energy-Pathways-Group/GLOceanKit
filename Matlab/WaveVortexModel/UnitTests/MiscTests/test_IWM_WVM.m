% IWVM decomposition test
%
% This code tests to make sure the decomposition by InternalWaveModel in 
% GLOceanKit gives the same results as the new WaveVortexModel for Winters model output
% 
% Output to command window:
% ~ IWM = max difference between u (from winter's output) and u_minus + u_plus + u_g (from InternalWaveModel)
% ~ WVM = max difference between u (from winter's output) and u_minus + u_plus + u_g (from WaveVortexModel)
% ~ diff_u_minus = max difference between u_minus (from InternalWaveModel) and u_minus (from WaveVortexModel)
% ~ diff_u_plus = max difference between u_plus (from InternalWaveModel) and u_plus (from WaveVortexModel)
% ~ diff_u_g = max difference between u_g (from InternalWaveModel) and u_g (from WaveVortexModel)
% ~ diff_u_full = max difference between u_plus + u_minus + u_g (from InternalWaveModel) and ...
%                                u_plus + u_minus + u_g (from WaveVortexModel)
%
% Created by Bailey Avila        12/02/2021
% Last modified                  01/14/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup!
clear all
addpath(genpath('~/matlab'));
addpath(genpath('../GLOceanKit-master'));
addpath(genpath('../GLNumericalModelingKit-master'));

movieflag = 0;
printflag = 0;

basedir = '/home/bavila/IWVM/';
%basedir = '/usr3/projects/IWVM/';

runroot = 'GM_500_01r07IWr01IW';

runDIR = [basedir,'model_raw/',runroot,'/'];

plotDIR = [basedir,'model_processed/',runroot,'/'];         % assign directory to save plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start by getting list of file names
fnms3D = dir([runDIR,'3D/*000.nc']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load problem params, io params, and high order diffusion params
pp = load_problem_params([runDIR,runroot,'_setup/']);
ip = load_io_params([basedir,'model_raw/',runroot,'/codes_etc']);
hodp = load_hod_params([basedir,'model_raw/',runroot,'/codes_etc']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the GLOceanKit wavemodel
% Jeffrey's code needs latitude, but it's not used in Winters code ... just calculate it from problem_params
latitude = asind(pp.f(1)/2/7.2921E-5);

% load in the EarlyIC_params - N0 set within here
run([runDIR,runroot,'_setup/input/EarlyIC_params']);

% initialize InternalWaveModel and WaveVortexModel
wavemodel = InternalWaveModelConstantStratification([pp.Lx,pp.Ly,pp.Lz],[pp.nx,pp.ny,pp.nz], latitude, N0, pp.rho_0);
wvm = WaveVortexModelConstantStratification([pp.Lx,pp.Ly,pp.Lz],[pp.nx,pp.ny,pp.nz], latitude, N0, pp.rho_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set some axis and colorbar limits
xlims = pp.Lx/2+pp.Ly/2*[-1 1];
ylims = [0 pp.Ly];
zlims = pp.Lz/2+pp.Lz/4*[-1 1];

% full 3D field dimensions
dims = [pp.Lx,pp.Ly,pp.Lz];

n = 501; % choosing a time step to look at from winters output
[x,y,z,time,u,v,w,s1,s2,phi,s1_bar,s2_bar] = load_XYZ_file(runroot,fnms3D(n).name(1:end-7)); %open/run load_XYZ_file

% set grid sizes and 3d X, Y, Z and other things only computed once
[Y, X, Z] = meshgrid(y, x, z);				  % for regular 3-D plots
dx = diff(x(1:2));                             % (m)
dy = diff(y(1:2));                             % (m)
dz = diff(z(1:2));                             % (m)

[~,~,S1_BAR] = meshgrid(ones(1,pp.ny),ones(1,pp.nx),s1_bar); %turning 1d into 3d
[dS1_BAR_dy,dS1_BAR_dx,dS1_BAR_dz] = gradient(S1_BAR,dy,dx,dz); % d/dz (rho_bar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InternalWaveModel decomposition
% decompose the fields into wave and vortical mode components
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(time, u, v, s1);

phase_plus = exp(sqrt(-1)*wavemodel.Omega*time);
phase_minus = exp(-sqrt(-1)*wavemodel.Omega*time);

u_wavep = wavemodel.transformToSpatialDomainWithF(wavemodel.u_plus.*phase_plus);	% wave plus portion of u
u_wavem = wavemodel.transformToSpatialDomainWithF(wavemodel.u_minus.*phase_minus);	% wave minus portion of u
u_g = wavemodel.u_g;								% vortical portion of u

v_wavep = wavemodel.transformToSpatialDomainWithF(wavemodel.v_plus.*phase_plus);	% wave plus portion of v
v_wavem = wavemodel.transformToSpatialDomainWithF(wavemodel.v_minus.*phase_minus);	% wave minus portion of v
v_g = wavemodel.v_g;								% vortical portion of v

w_wavep = wavemodel.transformToSpatialDomainWithG(wavemodel.w_plus.*phase_plus);	% wave plus portion of w
w_wavem = wavemodel.transformToSpatialDomainWithG(wavemodel.w_minus.*phase_minus);	% wave minus portion of w
w_g = 0*wavemodel.v_g;								% vortical portion of w=0 by definition

zeta_wavep = wavemodel.transformToSpatialDomainWithG(wavemodel.zeta_plus.*phase_plus);	% wave plus portion of zeta
zeta_wavem = wavemodel.transformToSpatialDomainWithG(wavemodel.zeta_minus.*phase_minus);	% wave minus portion of zeta
zeta_g = wavemodel.zeta_g;									% vortical portion of zeta

s1_wavep = zeta_wavep / (wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2));             % wave plus portion of s1
s1_wavem = zeta_wavem / (wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2));             % wave minus portion of s1
s1_g = wavemodel.zeta_g / (wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2));           % vortical portion of s1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WaveVortexModel decomposition
% decompose the fields into wave and vortical mode components using new WaveVortexModel
eta = s1 ./ dS1_BAR_dz;  % eta = rho'/(drhobar/dz) per Early et al 2021
[wvm.Ap,wvm.Am,wvm.A0] = wvm.transformUVEtaToWaveVortex(u,v,eta);

% zero out certain amplitudes from the WaveVortexModel decomposition to return either plus, minus, or geostrophic components of fields
[u_wavep_wv,v_wavep_wv,w_wavep_wv,zeta_wavep_wv] = wvm.transformWaveVortexToUVWEta(wvm.Ap,0,0); % plus component of fields
[u_wavem_wv,v_wavem_wv,w_wavem_wv,zeta_wavem_wv] = wvm.transformWaveVortexToUVWEta(0,wvm.Am,0); % minus component of fields
[u_g_wv,v_g_wv,w_g_wv,zeta_g_wv] = wvm.transformWaveVortexToUVWEta(0,0,wvm.A0); % geostrophic component of fields

s1_wavep = zeta_wavep_wv / (wvm.g/(wvm.rho0 * wvm.N0.^2));             % wave plus portion of s1
s1_wavem = zeta_wavem_wv / (wvm.g/(wvm.rho0 * wvm.N0.^2));             % wave minus portion of s1
s1_g = zeta_g_wv / (wvm.g/(wvm.rho0 * wvm.N0.^2));           % vortical portion of s1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots below are all pcolor plots of the u velocity field
% these are x-z plots looking at middle y location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting decomposed u fields using InteralWaveModel
figure(1)
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot u (before decomposition)
subplot(2,2,1)
pcolor(x/1000,z,squeeze(u(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u (m/s)')

colorbar
%caxis([-1*UAmp,1*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot (u_minus + u_plus) from InternalWaveModel
subplot(2,2,2)
pcolor(x/1000,z,squeeze(u_wavem(:,pp.ny/2,:)+u_wavep(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^+ + u^- (m/s)')

colorbar
%caxis([-0.6*UAmp,0.6*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot u_g from InternalWaveModel
subplot(2,2,4)
pcolor(x/1000,z,squeeze(u_g(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^0 (m/s)')

colorbar
%caxis([-0.5*UAmp,0.5*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot difference between u and (u_plus + u_minus + u_g) from InternalWaveModel
subplot(2,2,3)
pcolor(x/1000,z,squeeze(u(:,pp.ny/2,:)-(u_wavem(:,pp.ny/2,:)+u_wavep(:,pp.ny/2,:)+u_g(:,pp.ny/2,:)))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u-(u^++u^-+u^0) (m/s)')
IWM = max(abs((u(:)-(u_wavem(:)+u_wavep(:)+u_g(:)))))
colorbar
%caxis([-2e-5 2e-5])
bigtitle({'InternalWaveModel'})
if printflag
  print('-djpeg',[plotDIR,'IWVM_decomp_u.jpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting decomposed u fields using WaveVortexModel
figure(2)
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot u field (before decomposing)
subplot(2,2,1)
pcolor(x/1000,z,squeeze(u(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u (m/s)')

colorbar
%caxis([-1*UAmp,1*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot (u_minus + u_plus) from WaveVortexModel
subplot(2,2,2)
pcolor(x/1000,z,squeeze(u_wavem_wv(:,pp.ny/2,:)+u_wavep_wv(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^+ + u^- (m/s)')

colorbar
%caxis([-0.6*UAmp,0.6*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot u_g from WaveVortexModel
subplot(2,2,4)
pcolor(x/1000,z,squeeze(u_g_wv(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^0 (m/s)')

colorbar
%caxis([-0.5*UAmp,0.5*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot difference between u and (u_plus + u_minus + u_g) from WaveVortexModel
subplot(2,2,3)
pcolor(x/1000,z,squeeze(u(:,pp.ny/2,:)-(u_wavem_wv(:,pp.ny/2,:)+u_wavep_wv(:,pp.ny/2,:)+u_g_wv(:,pp.ny/2,:)))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u-(u^++u^-+u^0) (m/s)')
WVM = max(abs((u(:)-(u_wavem_wv(:)+u_wavep_wv(:)+u_g_wv(:)))))
colorbar
%caxis([-2e-5 2e-5])

bigtitle({'WaveVortexModel'})
if printflag
  print('-djpeg',[plotDIR,'IWVM_decomp_u_wv.jpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting difference between decomposed fields InternalWaveModel - WaveVortexModel
figure(3)
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot difference between u_plus from InternalWaveModel and from WaveVortexModel
subplot(2,2,1)
pcolor(x/1000,z,squeeze(u_wavep(:,pp.ny/2,:)-u_wavep_wv(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^+ (m/s)')
diff_u_plus = max(max(abs(u_wavep(:,pp.ny/2,:)-u_wavep_wv(:,pp.ny/2,:))))
colorbar
%caxis([-1*UAmp,1*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot difference between u_minus from InternalWaveModel and from WaveVortexModel
subplot(2,2,2)
pcolor(x/1000,z,squeeze(u_wavem(:,pp.ny/2,:) - u_wavem_wv(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^- (m/s)')
diff_u_minus = max(max(abs(u_wavem(:,pp.ny/2,:)-u_wavem_wv(:,pp.ny/2,:))))
colorbar
%caxis([-0.6*UAmp,0.6*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot difference between u_g from InternalWaveModel and from WaveVortexModel
subplot(2,2,4)
pcolor(x/1000,z,squeeze(u_g(:,pp.ny/2,:) - u_g_wv(:,pp.ny/2,:))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('u^0 (m/s)')
diff_u_g = max(max(abs(u_g(:,pp.ny/2,:)-u_g_wv(:,pp.ny/2,:))))
colorbar
%caxis([-0.5*UAmp,0.5*UAmp])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot difference between (u_plus + u_minus + u_g) from InternalWaveModel and from WaveVortexModel
subplot(2,2,3)
pcolor(x/1000,z,squeeze((u_wavem(:,pp.ny/2,:)+u_wavep(:,pp.ny/2,:)+u_g(:,pp.ny/2,:))...
  -(u_wavem_wv(:,pp.ny/2,:)+u_wavep_wv(:,pp.ny/2,:)+u_g_wv(:,pp.ny/2,:)))');
shading flat

xlabel('x (km)')
ylabel('z (m)')
title('(u^++u^-+u^0) (m/s)')
diff_u_full = max(max(abs((u_wavem(:,pp.ny/2,:)+u_wavep(:,pp.ny/2,:)+u_g(:,pp.ny/2,:))...
  -(u_wavem_wv(:,pp.ny/2,:)+u_wavep_wv(:,pp.ny/2,:)+u_g_wv(:,pp.ny/2,:)))))
colorbar
%caxis([-2e-5 2e-5])

bigtitle({'InternalWaveModel - WaveVortexModel'})
if printflag
  print('-djpeg',[plotDIR,'IWVM_decomp_u_diff.jpg']);
end
