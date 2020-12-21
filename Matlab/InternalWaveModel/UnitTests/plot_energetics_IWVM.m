function plot_energetics(runroot,varargin)
% function plot_energetics(runroot,varargin)
%
% INPUT:
%   runroot = 'GM_test_10';		-- optional (note, do not include "_setup" suffix, detects automatically)
%
% OPTIONAL INPUT:
%   'printflag',value   - 1=print figure to file, 0=do not
%   'saveflag',value    - 1=save diffusivity results to file, 0=do not
%   'fnmss',value	- increment with which to subsample .nc filenames
%
% Plots spectra and energy time series for BOUS model runs
%
% Written by Miles A. Sundermeyer, 1/31/18
% Last modified by Bailey Avila, 12/13/20 to add reference spectrum, and 
% change size of z direction (pp.nz-2 instead of pp.nz-1)

Version = 'plot_energetics_IWVM.m 12/4/2020';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for testing - skip function definition
if(~exist('runroot'))
  runroot = 'GM_500_01'
  varargin = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some flags - these values are overridden if passed as arguments to the function
printflag = 1;                          % print figures to file
saveflag = 1;                           % print figures to file
fnmss = 1;                              % subsample netcdf files by this amount

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% override any flags if they were passed as varargin
if length(varargin)==0
  % do nothing - go with above defaults
elseif length(varargin)>0
  nargs = length(varargin);
  if mod(nargs,2) ~= 0
    error('Arguments beyond runroot must be given as name/value pairs');
  end
  for nn = 1:2:nargs
    if ischar(varargin{nn+1})
      eval([varargin{nn},' = ''',varargin{nn+1},''';']);
    else
      eval([varargin{nn},' = ',num2str(varargin{nn+1}),';']);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up directories and get run info
basedir = '/usr3/projects/IWVM/';
basedir = '/home/bavila/IWVM/';
procDIR = [basedir,'model_processed/',runroot,'/'];
refDIR = [basedir,'model_raw/GM_500_01r06/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make model_processed directory if it does not already exist
if ~exist(procDIR)==7
  mkdir(procDIR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove and re-add my matlab directory, GLOceanKit, and GLOceanModelingKit in case anything has changed since last loaded
rmpath(genpath( [basedir 'GLOceanKit-master/Matlab'] ))
addpath(genpath( [basedir 'GLOceanKit-master/Matlab'] ))

rmpath(genpath( [basedir 'GLNumericalModelingKit-master/Matlab'] ))
addpath(genpath( [basedir 'GLNumericalModelingKit-master/Matlab'] ))

rmpath(genpath('~/matlab'))
addpath(genpath('~/matlab'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initiate some figures that will have cumulative info over all time steps
numfigs = 4;
if(0)					% open figures as we go
  for nn=1:numfigs
    figHandle(nn) = figure;
  end
  
  disp(' Please position figures now for optimal viewing ;) ...')
  mypause
else
  figHandle = [1:numfigs];		% force to use particular figure numbers
  for nn=1:numfigs
    figure(figHandle(nn))
    clf
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if this is a setup directory, or if we've already run the model
if exist([basedir,'model_raw/',runroot,'_setup/'])==7;	% still setup directory
  runDIR = [basedir,'model_raw/',runroot,'_setup/'];
  inputsub = 'input/';
  
  runningflag = input(' Setup directory still exists ... is the model still running? (yes=1,no=0)[1]: ');
  if isempty(runningflag)
    runningflag = 1;
  end
  
  gzflag = 0;				% save whether to gzip files again
  
  try                                   % presume IC filename and get rho_target, if it exists
    icfnm = [runDIR,inputsub,'SaveIC_EarlyIWmodel.nc'];
    rho_target = ncreadatt(icfnm,'/','rho_target');
  catch
  end
  
elseif exist([basedir,'model_raw/',runroot,'/'])==7;	% model already run, no "_setup" suffix
  runDIR = [basedir,'model_raw/',runroot,'/'];
  inputsub = 'codes_etc/input/';
  
  runningflag = 0;
  gzflag = 1;				% save whether to gzip files again
  disp(' Unzipping files from input directory ... (this could take a while for .nc file) ...')
  eval(['!gunzip ',runDIR,inputsub,'problem_params']);	% unzip to open
  eval(['!gunzip ',runDIR,inputsub,'io_params']);	% unzip to open
  eval(['!gunzip ',runDIR,inputsub,'SaveIC_EarlyIWmodel.nc']);	% unzip to open
  disp(' Done unzipping files from input directory ... (this could take a while for .nc file) ...')
  
  try                                   % presume IC filename and get rho_target, if it exists
    icfnm = [runDIR,runroot,'_setup/input/SaveIC_EarlyIWmodel.nc'];
    rho_target = ncreadatt(icfnm,'/','rho_target');
  catch
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load problem_params file to go along with this run
pp = load_problem_params([runDIR,inputsub(1:end-6)]);
ip = load_io_params([runDIR,inputsub(1:end-6)]);
run([runDIR,inputsub(1:end-6),'input/EarlyIC_params'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through output files plotting energy and spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start by getting list of file names
if runningflag
  fnms3D = dir([basedir,'BOUS/output/3D/*000.nc']);
else
  fnms3D = dir([runDIR,'3D/*000.nc']);
end

if(fnmss>1)
  disp(' Subsampling movie frames ...')
  %mypause
  fnms3D = fnms3D(1:fnmss:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set some axis and colorbar limits
xlims = pp.Lx/2+pp.Ly/2*[-1 1];
ylims = [0 pp.Ly];
zlims = [0 pp.Lz];
zlims = pp.Lz/2+pp.Lz/4*[-1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set some colors for spectra and Ri# histograms
colorset = colormap(jet(length(fnms3D)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load reference data
if ~strcmp(runroot,'GM_500_01r06')
  refs = load('/home/bavila/IWVM/model_raw/GM_500_01r06/plot_energetics_IWVM_new');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear some arrays
time = nan*[1:length(fnms3D)];

% kinetic and potential energy time series
KE = time;
PE = time;
KE_realspace = time;
APE_realspace = time;

% wave-vortex decomposition time series
ewave = time;
evortex = time;
eshear = time;
kewave = time;
kevortex = time;
pewave = time;
pevortex = time;

kwavplot = ones(pp.nx/2+1,1);
EA_isotr = nan*ones(pp.nx/2+1,length(time));
EG_isotr = EA_isotr;
KE_isotr = EA_isotr;
PE_isotr = EA_isotr;

wavem = ones(pp.nz-1,1);
vortexKE_m = nan*ones(pp.nz-2,length(time));
vortexPE_m = vortexKE_m;
wavePE_m = vortexKE_m;
wavePE_m = vortexKE_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' Reading through output files ...']);
for n=1:length(fnms3D)			% MAIN LOOP THROUGH .nc FILES
%for n=1%length(fnms3D)	%length(fnms3D)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp([' ',num2str(n),' of ',num2str(length(fnms3D))]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % assume data in slabs, load the slabs from each file and stack
  if runningflag
    error('This now broken - need to allow for BOUS/output path in `load_XYZ_file`');
  else
    [x_old,y_old,z_old,time(n),u,v,w,s1,s2,phi,s1_bar,s2_bar] = load_XYZ_file(runroot,fnms3D(n).name(1:end-7));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dx = pp.Lx/pp.nx;
  dy = pp.Ly/pp.ny;
  dz = pp.Lz/(pp.nz-1);
  
  % reset x,y,z arrays from scratch to avoid finite precision issues
  x = dx*[0:pp.nx-1]';
  y = dy*[0:pp.ny-1]';
  z = dz*[0:pp.nz-1]';
  
  % set up a 3D grid
  [Y,X,Z] = meshgrid(y,x,z);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate d/dz(s1_bar) for computing APE
  ds1_bar_dz = -N0.^2*(pp.rho_0/pp.g);
  
  % and make 3D array of s1_bar (not needed for ds1_bar_dz, since it is a single value under constant stratification)
  s1_BAR = permute(repmat(s1_bar,1,pp.nx,pp.ny),[2,3,1]);
  
  % make sure s1 is actually just rho_prime, not full rho_bar + rho_prime
  if(ip.write_s1_bar(1)==2)                             % s1 = rho_prime + rho_bar, subtract so s1 = rho_prime only
    s1 = s1 - s1_BAR;
  elseif(ip.write_s1_bar(1)==1)                         % s1 = rho_prime
    % do nothing
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set up Jeffrey's IW model
  % Jeffrey's code needs latitude, but it's not used in Winters code ... just calculate it from problem_params
  latitude = asind(pp.f(1)/2/7.2921E-5);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize Jeffrey's IW model based on what we used for model run i.c.
  wavemodel = InternalWaveModelConstantStratification([pp.Lx,pp.Ly,pp.Lz],[pp.nx,pp.ny,pp.nz], latitude, N0, pp.rho_0);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % After initializing the model with the correct dimensions and stratification, simply call either of the following.
  % Decomposed fields, Amp_plus, Amp_minus, B, and B0 are then contained in the wavemodel class structure.
  disp(' ... decomposing ...')
  %model.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields( time(n), u, v, zeta)
  wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(time(n), u, v, s1);	% Note: s1 is rho_prime w/o rho_bar
  zeta = s1 * wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2);
  disp(' ... done decomposing ...')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % These all double checked - sum of decomposed fields equals original fields, see test_IWVM_decomp.m
  if(0)			% use IW/VM decomposed fields to check we're doing things right
    % get wave and vortex pieces of each variable
    phase_plus = exp(sqrt(-1)*wavemodel.Omega*time(n));
    phase_minus = exp(-sqrt(-1)*wavemodel.Omega*time(n));
    
    u_wavep = wavemodel.TransformToSpatialDomainWithF(wavemodel.u_plus.*phase_plus);	% wave plus portion of u
    u_wavem = wavemodel.TransformToSpatialDomainWithF(wavemodel.u_minus.*phase_minus);	% wave minus portion of u
    u_g = wavemodel.u_g;								% vortical portion of u
    
    v_wavep = wavemodel.TransformToSpatialDomainWithF(wavemodel.v_plus.*phase_plus);	% wave plus portion of v
    v_wavem = wavemodel.TransformToSpatialDomainWithF(wavemodel.v_minus.*phase_minus);	% wave minus portion of v
    v_g = wavemodel.v_g;								% vortical portion of v
    
    w_wavep = wavemodel.TransformToSpatialDomainWithG(wavemodel.w_plus.*phase_plus);	% wave plus portion of w
    w_wavem = wavemodel.TransformToSpatialDomainWithG(wavemodel.w_minus.*phase_minus);	% wave minus portion of w
    w_g = 0*wavemodel.v_g;								% vortical portion of w=0 by definition
    
    zeta_wavep = wavemodel.TransformToSpatialDomainWithG(wavemodel.zeta_plus.*phase_plus);	% wave plus portion of zeta
    zeta_wavem = wavemodel.TransformToSpatialDomainWithG(wavemodel.zeta_minus.*phase_minus);	% wave minus portion of zeta
    zeta_g = wavemodel.zeta_g;									% vortical portion of zeta
    
    s1_wavep = zeta_wavep / (wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2));             % wave plus portion of s1
    s1_wavem = zeta_wavem / (wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2));             % wave minus portion of s1
    s1_g = wavemodel.zeta_g / (wavemodel.g/(wavemodel.rho0 * wavemodel.N0.^2));           % vortical portion of s1
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do some energy sum tests
  if(0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%% JEFFREY'S VERSION %%%%%
    % See: GLOceanKit/Matlab/InternalWaveModel/UnitTests/WaveVortexDecompositionTest.m
    integratedEnergy = trapz(wavemodel.z,mean(mean( u.^2 + v.^2 + w.^2 + N0*N0*zeta.*zeta, 1 ),2 ) )/2;
    fprintf('total integrated energy: %f m^3/s\n', integratedEnergy);
    
    spectralEnergy = sum(sum(sum(wavemodel.Amp_plus.*conj(wavemodel.Amp_plus) + wavemodel.Amp_minus.*conj(wavemodel.Amp_minus))));
    fprintf('total spectral energy: %f m^3/s\n', spectralEnergy);
    % NOTE: THE ABOVE ARE NOT THE SAME BECAUSE THE 2ND CALC IS ONLY WAVE ENERGY, NOT WAVE PLUS VORTEX!!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P2_pm = wavemodel.Amp_plus.*conj(wavemodel.Amp_plus) + wavemodel.Amp_minus.*conj(wavemodel.Amp_minus);
    HKE = sum(sum(sum(  P2_pm .* (1 + pp.f(1)*pp.f(1)./(wavemodel.Omega.*wavemodel.Omega)) ...
      .* (N0*N0 - wavemodel.Omega.*wavemodel.Omega) / (2 * (N0*N0 - pp.f(1)*pp.f(1)) )  )));
    VKE = sum(sum(sum(  P2_pm .* (wavemodel.Omega.*wavemodel.Omega - pp.f(1)*pp.f(1)) / (2 * (N0*N0 - pp.f(1)*pp.f(1)) )  )));
    PE = sum(sum(sum(  P2_pm .* N0*N0 .* (wavemodel.Omega.*wavemodel.Omega - pp.f(1)*pp.f(1)) ...
      ./ (2 * (N0*N0 - pp.f(1)*pp.f(1)) * wavemodel.Omega.*wavemodel.Omega )  )));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('total spectral wave energy (HKE + VKE + PE) = E: (%f + %f + %f) = %f m^3/s\n', HKE, VKE, PE, HKE+VKE+PE);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Wavenumbers in Jeffrey's model are as follows:
  % wavemodel.k increases from 2*pi*[0:pp.nx/2-1]/pp.Lx, then from -2*pi*[1:pp.nx/2]/pp.Lx
  my_k = [ [2*pi*[0:1:pp.nx/2-1]/pp.Lx] [2*pi*[-pp.nx/2:1:-1]/pp.Lx] ]';
  my_l = [ [2*pi*[0:1:pp.ny/2-1]/pp.Ly] [2*pi*[-pp.ny/2:1:-1]/pp.Ly] ]';
  % wavemodel.j goes from 1:nz-1 (nz odd)
  my_j = [1:pp.nz-1]';
  % wavemodel.[K,L,J] are meshgrid versions of wavemodel.[k,l,j]
  % wavemodel.K2 = wavemodel.K.^2 + wavemodel.L.^2,  and wavemodel.Kh = sqrt(wavemodel.K2)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute various time series of energy first to be sure normalizations all make sense
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First compute real space KE and APE
  KE_realspace(n) = 0.5*sum(u(:).^2+v(:).^2+w(:).^2)/(pp.nx*pp.ny*pp.nz)*pp.rho_0;		% (m^2/s^2 * kg/m^3 = J/m^3)
  APE_realspace(n) = 0.5*pp.g*sum(s1(:).^2./abs(ds1_bar_dz))/(pp.nx*pp.ny*pp.nz);		% (m^2/s^2 * kg/m^3 = J/m^3)
  
  % APE from Kundu - must divide this by rho0 to get J/kg
  % Note: multiply/divide by rho_0 to go between J/m^3 and J/kg
  APE2_realspace(n) = 0.5*pp.g^2/pp.rho_0/N0^2 * sum(s1(:).^2)/(pp.nx*pp.ny*pp.nz)/pp.rho_0;	% (J/m^3 / kg/m^3 = J/kg)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now from wave/vortex decomposition
  % See Jeffrey's:
  % GLOceanKit/Matlab/InternalWaveModel/InternalWaveModel.m
  %   Amp_plus, Amp_minus	% amplitudes of wave +/- modes
  %   B 			% amplitude of geostrophic component (baroclinic, matches Amp_plus/minus)
  %   B0 			% amplitude of geostrophic component (barotropic, k_z=0)
  
  % See: GLOceanKit/Matlab/WintersModel/Examples/ShowWaveVortexEnergyVsTime.m
  %   Should get: wavemodel.Apm_HKE_factor + wavemodel.Apm_VKE_factor + wavemodel.Apm_PE_factor == 1.0
  %   Why not also: wavemodel.B_HKE_factor + wavemodel.B_PE_factor/(pp.nx*pp.ny*pp.nz) + wavemodel.B0_HKE_factor/pp.nz == 1.0
  Ap2 = wavemodel.Amp_plus .* conj(wavemodel.Amp_plus);						% includes HKE, VKE and PE
  Am2 = wavemodel.Amp_minus .* conj(wavemodel.Amp_minus);					% includes HKE, VKE and PE
  B2 = (wavemodel.B_HKE_factor + wavemodel.B_PE_factor) .* (wavemodel.B .* conj(wavemodel.B));	% no VKE in geostrophic mode
  B02 = (wavemodel.B0_HKE_factor).* (wavemodel.B0 .* conj(wavemodel.B0));		% no PE in barotropic geostrophic mode
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Jeffrey's units are depth-integrated energy (m^3/s^2) so divide by pp.Lz to make this depth averaged (m^2/s^2),
  % and then multiply by density to make (m^2/s^2 * kg/m^3 = J/m^3)
  
  ewave(n) = pp.rho_0/pp.Lz * sum(sum(sum((Ap2 + Am2),1),2),3);
  kewave(n) = pp.rho_0/pp.Lz * sum(sum(sum((wavemodel.Apm_HKE_factor + wavemodel.Apm_VKE_factor).*(Ap2 + Am2),1),2),3);
  pewave(n) = pp.rho_0/pp.Lz * sum(sum(sum(wavemodel.Apm_PE_factor .*(Ap2 + Am2 ),1),2),3);
 
  evortex(n) =  pp.rho_0/pp.Lz * sum(sum(sum(B2,3)+B02,2),1);
  kevortex(n) = pp.rho_0/pp.Lz * sum(sum(sum((wavemodel.B_HKE_factor ...
    .* (wavemodel.B .* conj(wavemodel.B)) + repmat(B02,1,1,pp.nz-2)),1),2),3);
  pevortex(n) = pp.rho_0/pp.Lz * sum(sum(sum((wavemodel.B_PE_factor ...
    .* (wavemodel.B .* conj(wavemodel.B))                          ),1),2),3);
  
  eshear(n) = nan;
  
  KE(n) = kewave(n) + kevortex(n);
  PE(n) = pewave(n) + pevortex(n);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compare my calc's w/ Jeffrey's versions - nearly machine precision agreement?
  if(0)
    WaveEnergy_plus(n) = sum(Ap2(:));
    WaveEnergy_minus(n) = sum(Am2(:));
    VortexEnergy(n) = sum(B2(:));
    Vortex0Energy(n) = sum(B02(:));
    
    evortex(n) - pp.rho_0/pp.Lz*(VortexEnergy(n) + Vortex0Energy(n))
    ewave(n) - pp.rho_0/pp.Lz*(WaveEnergy_plus(n) + WaveEnergy_minus(n))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute full arrays of wave and vortex KE - note, choosing here to do KE spectra only, not KE and PE spectra
  % See: GLOceanKit/Matlab/WintersModel/Examples/CheckHyperviscosityBoundaries.m		% HKE, VKE
  wavePE_full  = wavemodel.Apm_PE_factor  .* ( Ap2 + Am2 );		% (m^3/s^2) wave plus & minus PE
  waveHKE_full = wavemodel.Apm_HKE_factor .* ( Ap2 + Am2 );		% (m^3/s^2) wave plus & minus HKE
  waveVKE_full = wavemodel.Apm_VKE_factor .* ( Ap2 + Am2 );		% (m^3/s^2) wave plus & minus VKE
  waveKE_full = waveHKE_full + waveVKE_full;				% (m^3/s^2) wave plus & minus KE = HKE + VKE
  %keyboard
  % note: the following leaves out barotropic mode HKE, since mode 0 and has no representation in a vertical spectrum
  vortexPE_full  = wavemodel.B_PE_factor   * abs(wavemodel.B).^2;	% (m^3/s^2) vortex baroclinic PE (for vert spectra)
  vortexHKE_full = wavemodel.B_HKE_factor .* abs(wavemodel.B).^2;	% (m^3/s^2) vortex baroclinic HKE (for vert spectra)
  
  % this one includes barotopic mode HKE by adding it after we have vertically summed all baroclinic mode energy
  vortexHKE_2D = 1/pp.Lz * sum( (wavemodel.B0_HKE_factor .* abs(wavemodel.B).^2) ,3) + B02;	% (m^3/s^2) (for horiz spectra)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sum across all k and l wavenumbers to get energies as fcn of vertical mode
  % Jeffrey's units are depth-integrated energy (m^3/s^2) so divide by pp.Lz to make depth averaged (m^2/s^2),
  wavePE_m(:,n) = 1/pp.Lz * squeeze(sum(sum(wavePE_full,1),2));		% (m^2/s^2)
  waveKE_m(:,n) = 1/pp.Lz * squeeze(sum(sum(waveKE_full,1),2));		% (m^2/s^2)
  
  vortexPE_m(:,n) = 1/pp.Lz * squeeze(sum(sum(vortexPE_full,1),2));		% (m^2/s^2)
  vortexKE_m(:,n) = 1/pp.Lz * squeeze(sum(sum(vortexHKE_full,1),2));		% (m^2/s^2)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % last step for vertical wavenumber spectra is to turn vertical mode number into vertical wavenumber,
  % and collapse everything (i.e., integrate) in x,y
  %wavem = wavemodel.j * 2*pi/(2*pp.Lz);					% (rad/m)
  wavem = wavemodel.j/2 * 2*pi/pp.Lz;					% (rad/m)
  
  if(0)
    % verify Parseval's Theorem
    % DOESNT QUITE WORK YET FOR VERTICAL!!
    %disp(' Parseval`s Theorem: [sum((u(:).^2 + v(:).^2)) sum(ufx2(:)+vfx2(:))] ...')
    %disp([sum((u(:).^2 + v(:).^2))/2 (sum(waveKE_m(:)+vortexKE_m(:)))*pp.nx*pp.ny*pp.nz])
    disp([mean((u(:).^2 + v(:).^2))/2 (sum(waveKE_m(:)+vortexKE_m(:)))])
    % ... but close.  This says that waveKE_m and vortexKE_m have units m^2/s^2
    
    integratedEnergy_h = trapz(wavemodel.z,mean(mean( u.^2 + v.^2, 1 ),2 ) )/2;
    disp(integratedEnergy_h/pp.Lz)
    spectralEnergy_h = sum(sum(sum(wavemodel.Apm_HKE_factor.*(wavemodel.Amp_plus.*conj(wavemodel.Amp_plus) + wavemodel.Amp_minus.*conj(wavemodel.Amp_minus)))));
    disp(spectralEnergy_h/pp.Lz)
    
    disp([mean((u(:).^2 + v(:).^2))/2])
    disp(1/pp.Lz * sum(sum(sum(wavemodel.Apm_HKE_factor.*(Ap2 + Am2)))));
    disp(1/pp.Lz * sum(sum(sum(waveHKE_full))));
    disp(sum(waveKE_m));
    disp([(sum(waveKE_m(:)+vortexKE_m(:)))])
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sum in z to get energy across all vertical modes as function of horizontal wavenumber
  % Jeffrey's units are depth-integrated energy (m^3/s^2) so divide by pp.Lz to make depth averaged (m^2/s^2),
  wavePE_kl = 1/pp.Lz * squeeze(sum(wavePE_full,3));			% (m^2/s^2)
  waveKE_kl = 1/pp.Lz * squeeze(sum(waveKE_full,3));			% (m^2/s^2)
  
  vortexPE_kl = 1/pp.Lz * squeeze(sum(vortexPE_full,3));
  vortexKE_kl = 1/pp.Lz * (squeeze(sum(vortexHKE_full,3))+vortexHKE_2D);	% is this correct? - check with Jeffrey
  
  vortexKE_kl = 1/pp.Lz * (sum( (wavemodel.B_HKE_factor .* abs(wavemodel.B).^2), 3) + wavemodel.B0_HKE_factor .*  abs(wavemodel.B0).^2 );
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % and fold the horizontal ones per wavenumber symmetry
  wavePE_klfold = wavePE_kl(1:pp.nx/2+1,:);			% this has zero wavenumber to max pos plus the first negative
  wavePE_klfold(2:end-1,:) = 2*wavePE_klfold(2:end-1,:);	% double the energy in all but zero mode and nyquist
  wavePE_klfold = wavePE_klfold(:,1:pp.ny/2+1);			% and same in l wavenumber direction
  wavePE_klfold(:,2:end-1) = 2*wavePE_klfold(:,2:end-1);
  
  waveKE_klfold = waveKE_kl(1:pp.nx/2+1,:);			% this has zero wavenumber to max pos plus the first negative
  waveKE_klfold(2:end-1,:) = 2*waveKE_klfold(2:end-1,:);	% double the energy in all but zero mode and nyquist
  waveKE_klfold = waveKE_klfold(:,1:pp.ny/2+1);			% and same in l wavenumber direction
  waveKE_klfold(:,2:end-1) = 2*waveKE_klfold(:,2:end-1);
  
  vortexPE_klfold = vortexPE_kl(1:pp.nx/2+1,:);			% this has zero wavenumber to max pos plus the first negative
  vortexPE_klfold(2:end-1,:) = 2*vortexPE_klfold(2:end-1,:);	% double the energy in all but zero mode and nyquist
  vortexPE_klfold = vortexPE_klfold(:,1:pp.ny/2+1);		% and same in l wavenumber direction
  vortexPE_klfold(:,2:end-1) = 2*vortexPE_klfold(:,2:end-1);
  
  vortexKE_klfold = vortexKE_kl(1:pp.nx/2+1,:);			% this has zero wavenumber to max pos plus the first negative
  vortexKE_klfold(2:end-1,:) = 2*vortexKE_klfold(2:end-1,:);	% double the energy in all but zero mode and nyquist
  vortexKE_klfold = vortexKE_klfold(:,1:pp.ny/2+1);		% and same in l wavenumber direction
  vortexKE_klfold(:,2:end-1) = 2*vortexKE_klfold(:,2:end-1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % and set the horizontal wavenumber arrays
  kwavplot = abs(wavemodel.k(1:pp.nx/2+1));
  lwavplot = abs(wavemodel.l(1:pp.ny/2+1));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % last step for horizontal wavenumber spectra is to make the horizontal energy arrays into isotropic wavenumber spectra
  % will use k-wavenumber as my isotropic wavenumber
  dk = 2*pi/pp.Lx;
  for nn = 1:length(kwavplot);
    % get indices to wavenumbers within this dk
    ind = find( abs(wavemodel.Kh(:,:,1)) >= (kwavplot(nn)-dk/2) ...
      & abs(wavemodel.Kh(:,:,1)) < (kwavplot(nn)+dk/2) );
    EA_isotr(nn,n) = sum(wavePE_kl(ind) + waveKE_kl(ind));
    EG_isotr(nn,n) = sum(vortexPE_kl(ind) + vortexKE_kl(ind));
    KE_isotr(nn,n) = sum(vortexKE_kl(ind)) + sum(waveKE_kl(ind));
    PE_isotr(nn,n) = sum(vortexPE_kl(ind)) + sum(wavePE_kl(ind));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot these to cross check
  if(0)
    figure(99)
    clf
    subplot(2,2,1)
    loglog(kwavplot/(2*pi),2*pi./kwavplot.*sum(wavePE_klfold+waveKE_klfold,2),'b-')
    hold on
    loglog(lwavplot/(2*pi),2*pi./lwavplot.*sum(wavePE_klfold+waveKE_klfold,1)','r-')
    loglog(kwavplot/(2*pi),2*pi./kwavplot.*EA_isotr,'k-')
    grid on
    xlim([1e-5 1e-2])
    ylim([1e-7 1e2])
    xlabel('cycles/m')
    ylabel('[m^2 s^{-2}] / [cycle/m]')
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % THIS REMAINDER OF CODE BASED ON OLD PLOT_ENERGETICS CODE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now compute and plot horizontal wavenumber spectra - follows spec_170420.m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute wavenumbers, periodic case
  dkx = 1./pp.Lx;				% (1/m)
  kx = dkx*(1:pp.nx/2+1)';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute fft
  ufx_avg = zeros(pp.nx,1);
  vfx_avg = zeros(pp.nx,1);
  for kk=1:pp.nz
    for jj=1:pp.ny
      % fft along this line
      ufx = squeeze(fft(u(:,jj,kk)));
      vfx = squeeze(fft(v(:,jj,kk)));
      
      % add to average
      ufx_avg = ufx_avg + abs(ufx).^2/pp.nx;
      vfx_avg = vfx_avg + abs(vfx).^2/pp.nx;
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compact this using symmetry - index (1) is DC mode, (pp.nx/2+1) is Nyquist
  ufx_avg(2:pp.nx/2) = ufx_avg(2:pp.nx/2) + flipud(ufx_avg(pp.nx/2+2:pp.nx));
  vfx_avg(2:pp.nx/2) = vfx_avg(2:pp.nx/2) + flipud(vfx_avg(pp.nx/2+2:pp.nx));
  ufx_avg(pp.nx/2+2:pp.nx) = [];
  vfx_avg(pp.nx/2+2:pp.nx) = [];
  
  ufx2 = ufx_avg;
  vfx2 = vfx_avg;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % %%%%% MAS VERSION %%%%%
  % verify Parseval's Theorem
  %disp(' Parseval`s Theorem: [sum((u(:).^2 + v(:).^2)) sum(ufx2(:)+vfx2(:))] ...')
  %disp([sum((u(:).^2 + v(:).^2)) sum(ufx2(:)+vfx2(:))])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(figHandle(1))
  % don't clear - need successive spectra to over plot
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,1)
  loglog(kx,ufx2/pp.nx, 'color',colorset(n,:),'LineWidth',1,'Markersize',3)
  hold on
  title('|U(k)|^2')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2)
  loglog(kx,vfx2/pp.nx, 'color',colorset(n,:),'LineWidth',1,'Markersize',3)
  hold on
  title('|V(k)|^2')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(n==1)
    for nn=1:2
      subplot(2,2,nn)
      %axis tight
      xlim([1e-5 1e-2]);
      %ylim([1e-4 1e5]);
      ylim([1e-7 1e2]);
      grid on
      
      xlabel('k [cycles m^{-1}]')
      ylabel('energy [m^3 s^{-2}]')
      ylabel('energy [m^2 s^{-2}]/[cycles/m]')
      %legend('|U(k)|^2','|V(k)|^2')
    end
  end
  
  bigtitle([verbatim(runDIR)])
  
  
  if(1)
    disp(' Computing PV ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute relative and potential vorticity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % See: GLOceanKit_080619_MAS/Matlab/InternalWaveModel/UnitTests/ComputeErtelPV.m
    % compute 3 components of relative vorticity
    zeta_x = DiffFourier(y,w,1,2) - DiffCosine (z,v,1,3);	% w_y - v_z	% (1/s)
    zeta_y = DiffCosine (z,u,1,3) - DiffFourier(x,w,1,1);	% u_z - w_x	% (1/s)
    zeta_z = DiffFourier(x,v,1,1) - DiffFourier(y,u,1,2);	% v_x - u_y	% (1/s)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scale b so that it is meters (isopycnal height), with a sign difference
    % Note: per Jeffrey's model, s1 should be rho_prime, not including rho_bar
    b = -(wavemodel.g/wavemodel.rho0)*s1/N0/N0;				% (m) i.e., b = rho_prime/(drho_bar/dz)
    
    % derivatives of buoyancy are unitless
    b_x = DiffFourier(x,b,1,1);						% (unitless)
    b_y = DiffFourier(y,b,1,2);						% (unitless)
    b_z = DiffSine   (z,b,1,3);						% (unitless)
    
    % this should be 1, if we've done this correctly
    bbar_z = wavemodel.N2AtDepth(wavemodel.Z)/N0/N0;			% (unitless)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute components of full PV (normalized by N0^2)
    PV_x = zeta_x        .* b_x;						% (1/s)
    PV_y = zeta_y        .* b_y;						% (1/s)
    PV_z = (zeta_z + pp.f(1)) .* b_z  +  zeta_z .* bbar_z;		% (1/s)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute linear PV (normalized by N0^2)
    PV_linear = zeta_z.*bbar_z  +  pp.f(1)*b_z;				% (1/s)
    
    PV = PV_x + PV_y + PV_z;						% (1/s)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute PV my way - NOT SURE IF THIS IS CORRECT - DOES NOT AGREE W/ JEFFREY'S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho_x = DiffFourier(x,s1,1,1);					% (unitless)
    rho_y = DiffFourier(y,s1,1,2);					% (unitless)
    rho_z = DiffSine   (z,s1,1,3) - pp.rho_0/pp.g*N0^2;			% (unitless)
    
    %RelVort = (dvdx - dudy);
    % vertical component of PV from Benoit Cushman Roisin
    %PV = RelVort - drhodz * pp.f(1)*pp.g/N0^2/pp.rho_0;			% MAS - replaced this on 1/3/2020
    
    % Linear PV
    PVlinear_my = ( (zeta_z + pp.f(1)).*rho_z ) / pp.rho_0;
    
    % Ertel PV
    PV_my = ( zeta_x.*rho_x + zeta_y.*rho_y + (zeta_z + pp.f(1)).*rho_z ) / pp.rho_0;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following bit derived from reviewer comments to Brunner-Suzuki et al (2012),
    % to check wheter linear wave-vortex decomposition is valid
    % - check if micro vertical Froude # << 1
    % - check if macro vertical Froude # << 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the micro (or vorticity-based) vertical Froude number (e.g., see Waite and Bartello, 2006)
    omegax = (dwdy - dvdz);
    omegay = -(dwdx - dudz);
    Frz = sqrt((omegax.^2 + omegay.^2)/2) / N0;
    Frz_micro = mean(Frz(:));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(figHandle(2))
  clf			% make anew each time with latest length time series
  
  keep = find(~isnan(time));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,1,1)
  if(1)
    plot(time(keep)/86400,KE(keep),'b-')
    hold on
    plot(time(keep)/86400,PE(keep),'r-')
    plot(time(keep)/86400,KE(keep)+PE(keep),'k-')
  else
    plot(time(keep)/86400,kewave(keep)+kevortex(keep),'m-')
    hold on
    plot(time(keep)/86400,pewave(keep)+pevortex(keep),'m-')
    plot(time(keep)/86400,ewave(keep)+evortex(keep),'k--')
  end
  
  % overplot my real-space versions of the same energy time series
  if(0)
    plot(time(keep)/86400,KE_realspace(keep),'c-')
    hold on
    plot(time(keep)/86400,APE_realspace(keep),'m-')
    plot(time(keep)/86400,KE_realspace(keep)+APE_realspace(keep),'k--')
  end
  
  % write at end of each line what's plotted
  text((time(n)+25)/86400, KE_realspace(n),'KE')
  text((time(n)+25)/86400, APE_realspace(n),'PE')
  text((time(n)+25)/86400, KE_realspace(n)+APE_realspace(n),'E_{tot}')
  
  xlabel('time (days)')
  ylabel('E (J/m^3)')
  
  legend('KE','PE','Total','location','southwest')
  %legend boxoff
  %xlim([0 max(time)/86400])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,1,2)
  plot(time(keep)/86400, ewave(keep),'b-')
  hold on
  plot(time(keep)/86400, evortex(keep),'r-')
  %plot(time(keep)/86400, eshear(keep),'m-')
  %plot(time(keep)/86400, ewave(keep) + eshear(keep) + evortex(keep),'k-')
  plot(time(keep)/86400, ewave(keep) + evortex(keep),'k-')
  
  % write at end of each line what's plotted
  text((time(n)+25)/86400, ewave(n)*0.95,'wave')
  text((time(n)+25)/86400, evortex(n)*1.2,'vortex')
  %text((time(n)+25)/86400, eshear(n)*1.1,'shear')
  %text((time(n)+25)/86400, (ewave(n) + eshear(n) + evortex(n))*1.1,'E_{tot}')
  text((time(n)+25)/86400, (ewave(n) + evortex(n))*1.1,'E_{tot}')
  
  xlabel('time (days)')
  ylabel('E (J m^{-3})')
  
  %legend('Wave','Vortex','Shear','Total','location','southwest')
  legend('Wave','Vortex','Total','location','southwest')
  %legend boxoff
  %xlim([0 max(time)/86400])
  
  bigtitle([verbatim(runDIR)])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot plan view and vertical slices of PV, velocity, rho, dye
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(figHandle(3));
  clf			% make anew each time with latest pcolor images
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plan view
  
  % do dye only on isopycnal
  if ~exist('rho_target')
    rho_target = mean(mean(s1(:,:,(pp.nz-1)/2)));         % just use mean density across domain for this
  end
  
  [C_interp,z_interp,rho_interp] = get_isopycnal_C(x,y,z,s1+s1_BAR,s2,rho_target,'min');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,1)
  %imagesc(x/1000,y/1000,squeeze(s2(:,:,(pp.nz-1)/2))')		% mid depth slice
  imagesc(x/1000,y/1000,C_interp')				% target isopycnal
  
  shading flat
  colormap(jet)
  %colorbar
  
  %xlabel('x (km)')
  ylabel('y (km)')
  title('dye (ppb)')
  
  axis equal
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(ylims/1000)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,2)
  Vel_h = sqrt(u.^2+v.^2);
  imagesc(x/1000,y/1000,squeeze( sqrt(u(:,:,(pp.nz-1)/2).^2 + v(:,:,(pp.nz-1)/2).^2) )')
  
  shading flat
  colormap(jet)
  %colorbar
  
  %xlabel('x (km)')
  %ylabel('y (km)')
  set(gca,'ytick','')
  title('|u_h| (m/s)')
  
  axis equal
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(ylims/1000)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,3)
  imagesc(x/1000,y/1000,squeeze(s1(:,:,(pp.nz-1)/2))')
  
  shading flat
  colormap(jet)
  %colorbar
  
  %xlabel('x (km)')
  %ylabel('y (km)')
  set(gca,'ytick','')
  title('\rho^\prime (kg/m^3)')
  
  axis equal
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(ylims/1000)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,4)
  imagesc(x/1000,y/1000,squeeze(PV(:,:,(pp.nz-1)/2))')
  
  shading flat
  colormap(jet)
  %colorbar
  
  %xlabel('x (km)')
  %ylabel('y (km)')
  set(gca,'ytick','')
  title('Ertel PV (s^{-1})')
  
  axis equal
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(ylims/1000)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % vertical sections
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,5)
  imagesc(x/1000,z,squeeze(s2(:,pp.ny/2,:))')
  %imagesc(x/1000,z,squeeze(s2(:,pp.ny/2,:))'-squeeze(s1(1,1,:))*ones(size(x'))+pp.rho_0)	% subtract a linear background
  %imagesc(x/1000,z,squeeze(s2(:,pp.ny/2,:))'+cos(2*pi*z/(pp.Lz*2))*ones(size(x')))		% subtract a cosine background
  
  shading flat
  colormap(jet)
  %colorbar
  
  xlabel('x (km)')
  ylabel('z (m)')
  %title('dye (ppb)')
  
  axis square
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(zlims)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,6)
  Vel_h = sqrt(u.^2+v.^2);
  imagesc(x/1000,z,squeeze( sqrt(u(:,pp.ny/2,:).^2 + w(:,pp.ny/2,:).^2) )')
  
  shading flat
  colormap(jet)
  %colorbar
  
  xlabel('x (km)')
  %ylabel('y (km)')
  set(gca,'ytick','')
  %title('|u_h| (m/s)')
  
  axis square
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(zlims)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,7)
  imagesc(x/1000,z,squeeze(s1(:,pp.ny/2,:))')
  
  shading flat
  colormap(jet)
  %colorbar
  
  xlabel('x (km)')
  %ylabel('y (km)')
  set(gca,'ytick','')
  %title('\rho^\prime (kg/m^3)')
  
  axis square
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(zlims)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,4,8)
  imagesc(x/1000,z,squeeze(PV(:,pp.ny/2,:))')
  
  shading flat
  colormap(jet)
  %colorbar
  
  xlabel('x (km)')
  %ylabel('y (km)')
  set(gca,'ytick','')
  %title('Ertel PV (s^{-1})')
  
  axis square
  set(gca,'ydir','normal')
  xlim(xlims/1000)
  ylim(zlims)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %packboth(2,4)
  
  bigtitle([verbatim(runDIR)])
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot various spectra - have the following from above
  
  m3s2_flag = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(figHandle(4));
  % don't clear - need successive spectra to over plot
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,1)
  %loglog(kwavplot/(2*pi),squeeze(mean(EA_isotr,2)),'Color',colorset(n,:),'LineWidth',1)
  if(m3s2_flag)
    loglog(kwavplot/(2*pi),2*pi./kwavplot.*EA_isotr(:,n),'Color',colorset(n,:),'LineWidth',1)
    ylabel('Horizontal [m^3 s^{-2}]')
    ylabel('Horizontal [m^2 s^{-2}] / [cycle/m]')
  else
    loglog(kwavplot/(2*pi),EA_isotr(:,n),'Color',colorset(n,:),'LineWidth',1)
    ylabel('Horizontal [m^2 s^{-2}]/[cycle/m]')
  end
  hold on
  
  colormap(jet)
  
  %xlabel('k (cycles/m)')
  title('Wave Energy')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2)
  %loglog(kwavplot/(2*pi),squeeze(mean(EG_isotr,2)),'Color',colorset(n,:),'LineWidth',1)
  if(m3s2_flag)
    loglog(kwavplot/(2*pi),2*pi./kwavplot.*EG_isotr(:,n),'Color',colorset(n,:),'LineWidth',1)
  else
    loglog(kwavplot/(2*pi),EG_isotr(:,n),'Color',colorset(n,:),'LineWidth',1)
  end
  hold on
  
  colormap(jet)
  
  %xlabel('k (cycles/m)')
  title('Vortex Energy')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,3)
  %loglog(wavem/(2*pi),squeeze(mean(EA_isotr,1)),'Color',colorset(n,:),'LineWidth',1)
  if(m3s2_flag)
    loglog(wavem/(2*pi),2*pi./wavem.*(wavePE_m(:,n)+waveKE_m(:,n)),'Color',colorset(n,:),'LineWidth',1)
    ylabel('Vertical [m^3 s^{-2}]')
    ylabel('Vertical [m^2 s^{-2}]/[cycle/m]')
  else
    loglog(wavem/(2*pi),(wavePE_m+waveKE_m),'Color',colorset(n,:),'LineWidth',1)
    ylabel('Vertical [m^2 s^{-2}]')
  end
  hold on
  
  colormap(jet)
  
  %xlabel('m (cycles/m)')
  xlabel('Wavenumber [cycles/m]')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,4)
  %loglog(wavem/(2*pi),squeeze(mean(EG_isotr,1)),'Color',colorset(n,:),'LineWidth',1)
  if(m3s2_flag)
    loglog(wavem/(2*pi),2*pi./wavem.*(vortexPE_m(:,n)+vortexKE_m(:,n)),'Color',colorset(n,:),'LineWidth',1)
  else
    loglog(wavem/(2*pi),(vortexPE_m+vortexKE_m),'Color',colorset(n,:),'LineWidth',1)
  end
  hold on
  
  colormap(jet)
  
  %xlabel('m (cycles/m)')
  xlabel('Wavenumber [cycles/m]')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % overplot reference spectrum if on last timestep
  if (n==length(fnms3D) && ~strcmp(runroot,'GM_500_01r06'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,1)
    if(m3s2_flag)
      loglog(refs.kwavplot/(2*pi),2*pi./refs.kwavplot.*refs.EA_isotr(:,1),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
      ylabel('Horizontal [m^3 s^{-2}]')
      ylabel('Horizontal [m^2 s^{-2}] / [cycle/m]')
    else
      loglog(refs.kwavplot/(2*pi),refs.EA_isotr(:,1),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
      ylabel('Horizontal [m^2 s^{-2}]/[cycle/m]')
    end
    hold on
    
    colormap(jet)
    
    title('Wave Energy')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,2)
    if(m3s2_flag)
      loglog(refs.kwavplot/(2*pi),2*pi./refs.kwavplot.*refs.EG_isotr(:,1),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
    else
      loglog(refs.kwavplot/(2*pi),refs.EG_isotr(:,1),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
    end
    hold on
    
    colormap(jet)
    
    title('Vortex Energy')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,3)
    if(m3s2_flag)
      loglog(refs.wavem/(2*pi),2*pi./refs.wavem.*(refs.wavePE_m(:,1)+refs.waveKE_m(:,1)),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
      ylabel('Vertical [m^3 s^{-2}]')
      ylabel('Vertical [m^2 s^{-2}]/[cycle/m]')
    else
      loglog(refs.wavem/(2*pi),(refs.wavePE_m(:,1)+refs.waveKE_m(:,1)),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
      ylabel('Vertical [m^2 s^{-2}]')
    end
    hold on
    
    colormap(jet)
    
    xlabel('Wavenumber [cycles/m]')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,4)
    if(m3s2_flag)
      loglog(refs.wavem/(2*pi),2*pi./refs.wavem.*(refs.vortexPE_m(:,1)+refs.vortexKE_m(:,1)),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
    else
      loglog(refs.wavem/(2*pi),(refs.vortexPE_m+refs.vortexKE_m(:,1)),'Color',colorset(1,:),'LineWidth',1,'LineStyle','--')
    end
    hold on
    
    colormap(jet)
    
    xlabel('Wavenumber [cycles/m]')
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % adjust axes
  for nn=1:4
    subplot(2,2,nn)
    
    if (nn==1 | nn==2)
      xlim([1e-5 1e-2])
    else
      xlim([0.5e-3 1e-0])
    end
    
    ylim([1e-7 1e2])
    
    if strcmp(runroot,'GM_500_00')
      ylim([1e-15 1e5])
    end
    grid on
  end
  
  warning(' CHECK NORMALIZATIONS ON WAVE-VORTEX SPECTRA ')
  
  if n==length(fnms3D)
    %    keyboard
  end
  %  drawnow
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end				% MAIN LOOP THROUGH .nc FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printflag
  for nn=[1:numfigs]
    figure(figHandle(nn))
    signature('',verbatim('/usr3/projects/IWVM/analysis/plot_energetics_IWVM_new_'),[.05 0.02 .11 .11]);
    print('-djpeg',[runDIR,'plot_energetics_IWVM_new_',num2str(nn),'.jpg']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust axes of 4th figure to zoom out y axis
ylim([1e-14 1e2])
figure(figHandle(4))

for nn=1:3
  subplot(2,2,nn)
  ylim([1e-14 1e2])
 
  grid on
end

if printflag
  figure(figHandle(4))
  signature('',verbatim(['/home/bavila/IWVM/analysis/plot_energetics_IWVM_new_',runroot]),[.05 0.02 .11 .11]);
  print('-djpeg',[runDIR,'plot_energetics_IWVM_zoomout_',num2str(4),'.jpg']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveflag
  save([runDIR,'plot_energetics_IWVM_new'],'time','KE','PE','KE_realspace','APE_realspace','ewave','kewave','pewave','evortex','kevortex','pevortex','kwavplot','EG_isotr','EA_isotr','KE_isotr','PE_isotr','wavem','wavePE_m','waveKE_m','vortexPE_m','vortexKE_m','colorset')
end
