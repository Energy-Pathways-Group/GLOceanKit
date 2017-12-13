%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SaveIC_EarlyIWmodel.m
%
% This script generates an internal wave field using Jeffrey Early's 
% InternalWaveModel.m code. Initial conditions correspond to a random GM
% wave field in a domain of the specified size.  The script saves a netcdf
% file that can be read by flow_solve_user.f90 and used by the Winters
% internal wave model.
%
% Also saves vertical structures, frequency, and horizontal wavenumber for
% periodic forcing.  Forcing is implemented in flow_solve_user.f90.
%
% The density field is decomposed as rho(x,y,z,t) = rhobar(x,y,z) + rhoprime(x,y,z,t)
%
% INPUT
%   runDIR: directory of a model run, with problem_params file setup as desired
%   rhobar: time-independent density field (kg m^-3)
%           Also called s1_bar in flow_solve_user.f90
%   GMReferenceLevel: the relative GM amplitude
%   omega_forcing: vector of forcing frequencies
% 
% OUTPUT saved in netcdf file
%   u:      x velocity (m s^-1)
%   v:      y velocity (m s^-1)
%   w:      z velocity (m s^-1)
%   rhoprime:   density anomaly (kg m^-3)
%               also called s1 in flow_solve_user.f90
%   drhobar_dz: vertical derivative of rhobar
%   drhobar_dzdz: second vertical derivative of rhobar
%
% NOTE
%   -Still need to implement arbitrary stratification.
%   -Further, assuming rhobar is uniform over the domain.
%   -s2_bar and s2' are unused here, but would specify initial conditions
%   for a second scalar field.
%   -must specify path to InternalWaveModel, below
%   -It might be useful to add a note to the output.  Could use that to
%   describe the density field. (e.g. "exponential stratification with 500m
%   e-folding depth").
%
% C. Wortham, December 22, 2016
% March 24, 2017: Add forcing vertical structures to output. CW
% May 1, 2017: Now using load_problem_params.  CW
% May 16, 2017: Added Jeffery's particle position initialization. CW
% May 18, 2017: Added offset to saved z vector by pp.Lz to be consistent
% with Winters model.  Remember that winters model has z positive upward,
% with origin at the bottom. CW
% May 22: Corrected forcing amplitude so that it's based on forcing
% structure functions averaged over the full domain, not just in vertical.
% Also, added some checks for consistency with problem_params. CW
% June 1: Cleaned up some of the scaling checks. CW
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add path to necessary code (user)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set basic directory for research
% basedir = '/Users/pascale/Desktop/Generate_ICS/';
basedir = '/Users/jearly/Documents/ProjectRepositories/GLOceanKit/Matlab/StokesDrift/Scratch/';
% basedir = '/home/cwortham/research/';

% add GLOceanKit directory to path
% addpath(genpath( ['/Users/pascale/Dropbox/NSF_IWV/GLOceanKit/Matlab/'] ))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the run directory and read in problem params (user)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runDIR = basedir;
pp = load_problem_params(runDIR);

% Jeffrey's code needs latitude, but it's not used in Winters'
% problem_params.  Just calculate it from f.
latitude = asind(pp.f(1)/2/7.2921E-5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% if restarting, specify the restart file.  This is used for positioning
% floats.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restart = false; % boolean
% restart_file = [runDIR 'input/XYZ_034550_restart.nc']; % file path


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify mean density profile rhobar(z) (user)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating density from stratification.
% could just specify rhobar directly.
% N0 = (5.2e-3)/2;
% z = linspace(-pp.Lz,0,pp.nz);
% rhobar = -pp.rho_0/pp.g * cumtrapz(z, N0^2*ones(size(z)));
% rhobar = rhobar - rhobar(end) + pp.rho_0;

% for arbitrary stratification, have to specify a function handle.
z = linspace(-pp.Lz,0,pp.nz);
[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
rhoFunction = @(z) rhoFunc(z) - rhoFunc(max(zIn)) + pp.rho_0;
rhobar = rhoFunction(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the wave field (user)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use this for full GM spectrum:
% wave_type = 'GM';
% GMReferenceLevel = 1.0*pp.Lz/(2*1300);

% use this for single plane wave:
wave_type = 'plane';
k0 = 4; % k=0..Nx/2
l0 = 0; % l=0..Ny/2
j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
UAmp = 0.1089; % m/s
sign = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the forcing frequencies/amplitudes (user)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set frequencies at which to force, rad s^-1
omega_forcing = [pp.f(1), 2*pi/(12*3600)];

% set amplitude as relative weight of various modes.  Should sum to 1.  
% it's modified to actual amplitude below.
amp_forcing = 0*[0.5,0.5];

% set random phase offset for each forcing frequency
phase = 2*pi*rand(size(omega_forcing)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the forcing modes (internal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the vertical modes code
modes = InternalModes(rhoFunction,zIn,z,latitude);

% take advantage of the density derivatives offered
N2 = modes.N2;
drhobar_dz = modes.rho_z';
drhobar_dzdz = modes.rho_zz';

% compute the mode shape at each frequency
F = nan(length(z),length(omega_forcing));
G = F;
h = nan(size(omega_forcing));
for i = 1:length(omega_forcing)
    [F_temp,G_temp,h_temp] = modes.ModesAtFrequency(omega_forcing(i));
    % save lowest vertical mode F, G in array    
    F(:,i) = F_temp(:,1);
    G(:,i) = G_temp(:,1);
    h(i) = h_temp(1);
end

% compute forcing vertical structures, based on Early et al. manuscript
% draft eq. 29 and notebook for 3/23/2017.

% K, k, l from omega... from dispersion relation
K = sqrt( (omega_forcing.^2 - pp.f(1)^2)./(pp.g*h) );
l_over_k = 0; % ratio l/k.  Assume it's zero below.
k = sqrt( K.^2/(1+l_over_k^2) );
l = sqrt( K.^2 - k.^2 );

% time-independent part of forcing terms.  
% This simpler form applies only for l=0.  Would have to include
% extra terms for arbitrary wave direction.
u_forcing = bsxfun(@times, 1./sqrt(h), F);
v_forcing = bsxfun(@times, -pp.f(1)./omega_forcing./sqrt(h), F);
w_forcing = bsxfun(@times, K.*sqrt(h), G);
rho_forcing = bsxfun(@times, K.*sqrt(h)./omega_forcing, bsxfun(@times, drhobar_dz', G));

% full domain vectors
[X,Y,Z] = ndgrid( linspace(0,pp.Lx,pp.nx), linspace(0,pp.Ly,pp.ny), linspace(-pp.Lz,0,pp.nz));

% amplitude and energy partition between modes
% This is a simplified version of the Sugiyama et al. (2009) formulation
% with one amplitude factor (alpha=beta).  See notebook for 3/24/2017.
gamma = amp_forcing/sum(amp_forcing); % fraction of energy put into each mode
Gamma = 0.2; % Osborn (1980) dissipation constant
kappa_rho = 1e-5; % typical eddy diffusivity coefficient
% epsilon = kappa_rho * mean(N2)/Gamma; % energy dissipation rate
epsilon = 3.5e-10; % dissipation rate based on unforced run of cim_dev
for ii=1:length(omega_forcing)
    u_temp = bsxfun(@times, cos(k(ii)*X + l(ii)*Y), permute(u_forcing(:,ii), [3 2 1]));
    v_temp = bsxfun(@times, sin(k(ii)*X + l(ii)*Y), permute(v_forcing(:,ii), [3 2 1]));
    w_temp = bsxfun(@times, sin(k(ii)*X + l(ii)*Y), permute(w_forcing(:,ii), [3 2 1]));
    rho_temp = bsxfun(@times, cos(k(ii)*X + l(ii)*Y), permute(rho_forcing(:,ii), [3 2 1]));
    amp_forcing(ii) = gamma(ii)*pp.rho_0*epsilon./ ( pp.rho_0* (mean(mean(u_temp(:,:,ii).^2,1),2) + mean(mean(v_temp(:,:,ii).^2,1),2) + mean(mean(w_temp(:,:,ii).^2,1),2) ) + pp.g^2*mean(mean(rho_temp(:,:,ii).^2,1),2)/pp.rho_0/mean(N2) );
end

% check that amplitudes satisfy Sugiyama et al (2009) eq. 6
% for ii=1:length(omega_forcing)
%     KEI = amp_forcing(ii)*pp.rho_0*(mean(u_forcing(:,ii).^2) + mean(v_forcing(:,ii).^2) + mean(w_forcing(:,ii).^2));
%     PEI = amp_forcing(ii)*pp.g^2*mean(rho_forcing(:,ii).^2)/pp.rho_0/mean(N2);
%     KEI/PEI
%     (omega_forcing(ii)^2+f0^2)/(omega_forcing(ii)^2-pp.f(1)^2)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model (internal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the model
% wavemodel = InternalWaveModelConstantStratification([pp.Lx, pp.Ly, pp.Lz], [pp.nx, pp.ny, pp.nz], latitude, N0);

% create a wave realization
if strcmp( wave_type, 'GM')
    % GM spectrum
    wavemodel = InternalWaveModelArbitraryStratification([pp.Lx, pp.Ly, pp.Lz], [pp.nx, pp.ny, pp.nz], rhoFunction, z, pp.nz-10, latitude);
    wavemodel.InitializeWithGMSpectrum(GMReferenceLevel,0); % 1=shouldRandomizeAmplitude, 0 keeps same amplitude for all waves of given wavenumber
elseif strcmp( wave_type, 'plane')
    % single plane wave
    wavemodel = InternalWaveModelArbitraryStratification([pp.Lx, pp.Ly, pp.Lz], [pp.nx, pp.ny, pp.nz], rhoFunction, z, j0+1, latitude);
    period = wavemodel.InitializeWithPlaneWave(k0, l0, j0, UAmp, sign);
    
    if 1 % Add the next order correction term
        kk = wavemodel.k(k0+1);
        h = wavemodel.h(k0+1,l0+1,j0);
        omega = 2*pi/period;
        epsilon = UAmp/(omega/kk);
        
        zHR = linspace(min(z),max(z),5000)';
        im = InternalModesWKBSpectral(rhoFunction,[min(z) max(z)],zHR,latitude);
        im.normalization = Normalization.uMax;
        [F_hr_out,G_hr_out,~,~] = im.ModesAtWavenumber(kk);
        F_hr = F_hr_out(:,j0);
        G_hr = G_hr_out(:,j0);
        
        f = -im.rho_zz .* G_hr .* G_hr / wavemodel.internalModes.rho0;
        
        im.normalization = Normalization.kConstant;
        [~,G_2k_out,h_2k_out,~] = im.ModesAtWavenumber(2*kk);
        a = zeros(size(h_2k_out));
        for i=1:size(G_2k_out,2)
            a(i) = trapz(zHR,f.*G_2k_out(:,i));
        end
        gamma = ((h*h_2k_out./(h-h_2k_out))).*a;
        
        wavemodel.internalModes.normalization = Normalization.kConstant;
        wavemodel.internalModes.nModes = length(a);
        [F_2k_out,G_2k_out,~,~] = wavemodel.internalModes.ModesAtWavenumber(2*kk);
        
        wavemodel.internalModes.normalization = Normalization.uMax;
        [~,G_out,~,~] = wavemodel.internalModes.ModesAtWavenumber(kk);
        G = G_out(:,j0);
        
        X = wavemodel.X;
        Y = wavemodel.Y;
        Z = wavemodel.Z;
        RhoBarDz = -(wavemodel.rho0/9.81)*wavemodel.N2AtDepth(Z);
        
        Phi = repmat(permute((F_2k_out*(((h./h_2k_out).*gamma).')),[3 2 1]),size(X,1),size(X,2));
        Gamma = repmat(permute((G_2k_out*(gamma.')),[3 2 1]),size(X,1),size(X,2));
        Rho_zzTerm = repmat(permute((h*wavemodel.internalModes.rho_zz.*G.*G),[3 2 1]),size(X,1),size(X,2));
        
        u_correction = (UAmp * epsilon / 4) * cos(2*kk*X) .* Phi;
        w_correction = (UAmp*kk*h*epsilon/2) * sin(2*kk*X) .* Gamma;
        rho_correction = (h*epsilon*epsilon/4) * cos(2*kk*X) .* ( Rho_zzTerm +  RhoBarDz .* Gamma );
    end
    
else
    error('Invalid wave_type specified.  Must be "GM" or "plane".');
end

% run the model.
time = 0;
[u,v]=wavemodel.VelocityFieldAtTime(time);
[w,zeta] = wavemodel.VerticalFieldsAtTime(time);

% compute rhoprime since it isn't output from InternalWaveModel explicitly
DRHODBAR_DZ = repmat( reshape(drhobar_dz, [1,size(drhobar_dz)]) ,pp.nx,pp.ny);
rhoprime = - DRHODBAR_DZ .* zeta;

wavemodel.z = wavemodel.z+pp.Lz; % add Lz to z for compatibility with flow_solve z=0 level which is at the bottom.

if exist('u_correction','var')
   u = u + u_correction;
   w = w + w_correction;
   rhoprime = rhoprime + rho_correction;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters (user)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set number of particles
% N = 10; % Create an NxN grid of floats
% nLevels = 5; % Number of levels spanning [0 -pp.Lz/2]
% 
% % check that nparticles from problem_params matches specification here
% if ~( N^2*nLevels == pp.nparticles )
%     error('nFloats specified here must match nparticles specified in problem_params file.');
% end
% 
% % distribute particles evenlly. Must be on grid points.
% % dx = wavemodel.x(2)-wavemodel.x(1);
% % dy = wavemodel.y(2)-wavemodel.y(1);
% % x_float = (0:N-1)*dx;
% % y_float = (0:N-1)*dy;
% % z_float = (0:nLevels-1)*(-pp.Lz/(2*(nLevels-1)));
% dx = pp.Lx/4/(N-1); % 1/4 of domain
% dy = pp.Ly/4/(N-1); % 1/4 of domain
% dz = -pp.Lz/2/(nLevels-1); % 1/2 of domain
% x_float = (0:N-1)*dx;
% x_float = x_float - mean(x_float) + pp.Lx/2; % shift to center of domain
% y_float = (0:N-1)*dy;
% y_float = y_float - mean(y_float) + pp.Ly/2; % shift to center of domain
% z_float = (0:nLevels-1)*dz; % note that these are in Jeffrey's vertical coordinate (zero at surface)
% 
% % expand particle position arrays
% [x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
% x_float = reshape(x_float,[],1);
% y_float = reshape(y_float,[],1);
% z_float = reshape(z_float,[],1);
% nFloats = numel(x_float);
% 
% % Now let's place the floats along an isopycnal.
% interpolationMethod = 'spline';
% if restart
%     % use full density from the restart file
%     Rho = ncread(restart_file, 'rho'); 
% else
%     % use new density generated by wavemodel
%     Rho = bsxfun(@plus, rhoprime, permute(rhobar,[3,1,2]));
% end
% z_isopycnal = PlaceParticlesOnIsopycnal(x_float,y_float,z_float, X,Y,Z,Rho, drhobar_dz, interpolationMethod,1e-8);
% 
% 
% % Shift vertical position to flow_solve coordinates (zero at surface)
% % x_float, y_float, z_float contain float positions on isopycnals.
% z_float = z_isopycnal + pp.Lz;  

x_float = [0 0 0]+pp.Lx/2;
y_float = [0 0 0]+pp.Ly/2;
z_float = [-10; -250; -625]+5000;
nFloats = numel(x_float);

% check that nparticles from problem_params matches specification here
if ~( nFloats == pp.nparticles )
    error('nFloats specified here must match nparticles specified in problem_params file.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check consistency of ICs with scalings in problem_params (internal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Checking consistency of initial conditions with scalings specified in problem_params...')

% set maximum allowed ratio for scalings
max_ratio = 5;

% stratification
strat_ratio = median(abs(drhobar_dz)) / pp.dgrad;
% if strat_ratio<1
%     strat_ratio = 1/strat_ratio;
% end

% horizontal velocity
% uvel_ratio = sqrt(mean(u(:).^2+v(:).^2)) / pp.u0;
uvel_ratio = max(abs(u(:))) / pp.u0;
% if uvel_ratio<1
%     uvel_ratio = 1/uvel_ratio;
% end

% density anomaly
% dens_ratio = sqrt(mean(rhoprime(:).^2)) / pp.scalar_scale(1);
dens_ratio = max(abs(rhoprime(:))) / pp.scalar_scale(1);
% if dens_ratio<1
%     dens_ratio = 1/dens_ratio;
% end

% report
valid_scaling = true; % just keeping track of whether scalings are valid
if (strat_ratio>max_ratio) || (strat_ratio<1/max_ratio)
    warning('Ratio of median drhbobar_dz to problem_params stratification scale is %1.10f.  Proceeding, but make sure this is what you want.', strat_ratio)
    valid_scaling = false;
end
if (uvel_ratio>max_ratio) || (uvel_ratio<1/max_ratio)
%     warning('Ratio of RMS velocity to problem_params velocity scale is %1.1f.  Proceeding, but make sure this is what you want.', uvel_ratio)
    warning('Ratio of max velocity to problem_params velocity scale is %1.1f.  Proceeding, but make sure this is what you want.', uvel_ratio)
    valid_scaling = false;
end
if (dens_ratio>max_ratio) || (dens_ratio<1/max_ratio)
%     warning('Ratio of RMS rhoprime to problem_params s1 scale is %1.1f.  Proceeding, but make sure this is what you want.', dens_ratio)
    warning('Ratio of max rhoprime to problem_params s1 scale is %1.1f.  Proceeding, but make sure this is what you want.', dens_ratio)
    valid_scaling = false;
end
if valid_scaling
    disp('...initial conditions consistent with scalings specified in problem_params.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check if forcing wavelengths fit in domain (internal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Checking whether or not the forcing wavelengths fit in the model domain...')

% what fraction of a cycle can forcing be off by?
tolerance = .1; 

% how close do the zonal wavenumbers fit in the domain?
k_misfit = mod( pp.Lx./( 2*pi./k ), 1);

% how close do the meridional wavenumbers fit in the domain?
l_misfit = mod( pp.Ly./( 2*pi./l ), 1);

% report
no_misfit = true; % just keeping track of whether any wavenumbers don't fit.
for ii=1:length(omega_forcing)
    if k_misfit(ii)>tolerance && k_misfit(ii)<1-tolerance
        warning('Zonal wavenumber for mode %d does not fit well in the domain. \n I recommend making the domain an integer multiple of forcing wavelength. \n 2*pi/k=%1.0fm and Lx=%1.0fm', ii, 2*pi/k(ii), pp.Lx )
        no_misfit = false;
    end
end
for ii=1:length(omega_forcing)
    if l_misfit(ii)>tolerance && l_misfit(ii)<1-tolerance
        warning('Meridional wavenumber for mode %d does not fit well in the domain. \n I recommend making the domain an integer multiple of forcing wavelength. \n 2*pi/l=%1.0fm and Ly=%1.0fm', ii, 2*pi/l(ii), pp.Ly )
        no_misfit = false;
    end
end
if no_misfit
    disp('...all forcing wavenumbers fit well in the domain.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file (internal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename = sprintf('SaveIC_EarlyIWmodel_%dx%dx%d@%s.nc', pp.nx,pp.ny,pp.nz,datestr(datetime('now'),'yyyy-mm-dd'));
% filename = sprintf('SaveIC_EarlyIWmodel_%dx%dx%d.nc', pp.nx,pp.ny,pp.nz);
filename = 'SaveIC_EarlyIWmodel.nc';

ncid = netcdf.create(filename, 'CLOBBER');

precision = 'single';
if strcmp(precision,'single')
    ncPrecision = 'NC_FLOAT';
    setprecision = @(x) single(x);
    bytePerFloat = 4;
else
    ncPrecision = 'NC_DOUBLE';
    setprecision = @(x) double(x);
    bytePerFloat = 8;
end

% Define the dimensions
iDimID = netcdf.defDim(ncid, 'idimension', wavemodel.Nx);
jDimID = netcdf.defDim(ncid, 'jdimension', wavemodel.Ny);
kDimID = netcdf.defDim(ncid, 'kdimension', wavemodel.Nz);
tDimID = netcdf.defDim(ncid, 'timedimension', netcdf.getConstant('NC_UNLIMITED'));
fDimID = netcdf.defDim(ncid, 'forcedimension', length(omega_forcing));
pDimID = netcdf.defDim(ncid, 'particledimension', nFloats);

% Define the coordinate variables
iVarID = netcdf.defVar(ncid, 'x', ncPrecision, iDimID);
jVarID = netcdf.defVar(ncid, 'y', ncPrecision, jDimID);
kVarID = netcdf.defVar(ncid, 'z', ncPrecision, kDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
fVarID = netcdf.defVar(ncid, 'mode', ncPrecision, fDimID);
netcdf.putAtt(ncid,iVarID, 'long_name', 'x coordinate');
netcdf.putAtt(ncid,jVarID, 'long_name', 'y coordinate');
netcdf.putAtt(ncid,kVarID, 'long_name', 'z coordinate');
netcdf.putAtt(ncid,tVarID, 'long_name', 'time coordinate');
netcdf.putAtt(ncid,fVarID, 'long_name', 'forcing mode number');
netcdf.putAtt(ncid,iVarID, 'units', 'm');
netcdf.putAtt(ncid,jVarID, 'units', 'm');
netcdf.putAtt(ncid,kVarID, 'units', 'm');
netcdf.putAtt(ncid,tVarID, 'units', 's');
netcdf.putAtt(ncid,fVarID, 'units', 'unitless');

% Define the dynamical variables
uVarID = netcdf.defVar(ncid, 'u', ncPrecision, [iDimID,jDimID,kDimID,tDimID]);
vVarID = netcdf.defVar(ncid, 'v', ncPrecision, [iDimID,jDimID,kDimID,tDimID]);
wVarID = netcdf.defVar(ncid, 'w', ncPrecision, [iDimID,jDimID,kDimID,tDimID]);
rhoprimeVarID = netcdf.defVar(ncid, 's1', ncPrecision, [iDimID,jDimID,kDimID,tDimID]);
rhobarVarID = netcdf.defVar(ncid, 'rhobar', ncPrecision, kDimID);
drhobardzVarID = netcdf.defVar(ncid, 'drhobar_dz', ncPrecision, kDimID);
drhobardzdzVarID = netcdf.defVar(ncid, 'drhobar_dzdz', ncPrecision, kDimID);
% zetaVarID = netcdf.defVar(ncid, 'zeta', ncPrecision, [iDimID,jDimID,kDimID,tDimID]);
omegaForceVarID = netcdf.defVar(ncid, 'omega_forcing', ncPrecision, fDimID);
ampForceVarID = netcdf.defVar(ncid, 'amp_forcing', ncPrecision, fDimID);
phaseForceVarID = netcdf.defVar(ncid, 'phase_forcing', ncPrecision, fDimID);
kForceVarID = netcdf.defVar(ncid, 'k_forcing', ncPrecision, fDimID);
lForceVarID = netcdf.defVar(ncid, 'l_forcing', ncPrecision, fDimID);
uForceVarID = netcdf.defVar(ncid, 'u_forcing', ncPrecision, [kDimID, fDimID]);
vForceVarID = netcdf.defVar(ncid, 'v_forcing', ncPrecision, [kDimID, fDimID]);
wForceVarID = netcdf.defVar(ncid, 'w_forcing', ncPrecision, [kDimID, fDimID]);
rhoForceVarID = netcdf.defVar(ncid, 'rho_forcing', ncPrecision, [kDimID, fDimID]);
xfloatVarID = netcdf.defVar(ncid,  'x_float', ncPrecision, pDimID);
yfloatVarID = netcdf.defVar(ncid,  'y_float', ncPrecision, pDimID);
zfloatVarID = netcdf.defVar(ncid,  'z_float', ncPrecision, pDimID);
% assign long names
netcdf.putAtt(ncid,uVarID, 'long_name', 'u velocity');
netcdf.putAtt(ncid,vVarID, 'long_name', 'v velocity');
netcdf.putAtt(ncid,wVarID, 'long_name', 'w velocity');
netcdf.putAtt(ncid,rhoprimeVarID, 'long_name', 'time-dependent density anomaly');
netcdf.putAtt(ncid,rhobarVarID, 'long_name', 'time-independent density anomaly');
netcdf.putAtt(ncid,drhobardzVarID, 'long_name', 'time-independent density derivative');
netcdf.putAtt(ncid,drhobardzdzVarID, 'long_name', 'time-independent density second derivative');
% netcdf.putAtt(ncid,zetaVarID, 'units', 'm');
netcdf.putAtt(ncid,omegaForceVarID, 'long_name', 'forcing frequencies');
netcdf.putAtt(ncid,ampForceVarID, 'long_name', 'forcing amplitudes');
netcdf.putAtt(ncid,phaseForceVarID, 'long_name', 'forcing phase offset');
netcdf.putAtt(ncid,kForceVarID, 'long_name', 'forcing zonal wavenumber');
netcdf.putAtt(ncid,lForceVarID, 'long_name', 'forcing meridional wavenumber');
netcdf.putAtt(ncid,uForceVarID, 'long_name', 'vertical structure function for u');
netcdf.putAtt(ncid,vForceVarID, 'long_name', 'vertical structure function for v');
netcdf.putAtt(ncid,wForceVarID, 'long_name', 'vertical structure function for w');
netcdf.putAtt(ncid,rhoForceVarID, 'long_name', 'vertical structure function for rho');
netcdf.putAtt(ncid,xfloatVarID, 'long_name', 'float/particle zonal position');
netcdf.putAtt(ncid,yfloatVarID, 'long_name', 'float/particle meridional position');
netcdf.putAtt(ncid,zfloatVarID, 'long_name', 'float/particle vertical position');
% assign units
netcdf.putAtt(ncid,uVarID, 'units', 'm/s');
netcdf.putAtt(ncid,vVarID, 'units', 'm/s');
netcdf.putAtt(ncid,wVarID, 'units', 'm/s');
netcdf.putAtt(ncid,rhoprimeVarID, 'units', 'kg m^-3');
netcdf.putAtt(ncid,rhobarVarID, 'units', 'kg m^-3');
netcdf.putAtt(ncid,drhobardzVarID, 'units', 'kg m^-4');
netcdf.putAtt(ncid,drhobardzdzVarID, 'units', 'kg m^-5');
% netcdf.putAtt(ncid,zetaVarID, 'units', 'm');
netcdf.putAtt(ncid,omegaForceVarID, 'units', 'rad/s');
netcdf.putAtt(ncid,ampForceVarID, 'units', 's^-1');
netcdf.putAtt(ncid,phaseForceVarID, 'units', 'dimensionless');
netcdf.putAtt(ncid,kForceVarID, 'units', 'rad/m');
netcdf.putAtt(ncid,lForceVarID, 'units', 'rad/m');
netcdf.putAtt(ncid,uForceVarID, 'units', 'm/s');
netcdf.putAtt(ncid,vForceVarID, 'units', 'm/s');
netcdf.putAtt(ncid,wForceVarID, 'units', 'm/s');
netcdf.putAtt(ncid,rhoForceVarID, 'units', 'kg m^-3');
netcdf.putAtt(ncid,xfloatVarID, 'units', 'm');
netcdf.putAtt(ncid,yfloatVarID, 'units', 'm');
netcdf.putAtt(ncid,zfloatVarID, 'units', 'm');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'f0', pp.f(1));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'rho_0', pp.rho_0);
if strcmp(wave_type,'GM')
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
elseif strcmp(wave_type,'plane')
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'plane wave amplitude', UAmp);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'plane wave k0', k0);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'plane wave l0', l0);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'plane wave j0', j0);
else
    error('Invalid wave_type specified.  Must be "GM" or "plane".');
end
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'InternalWaveModel version', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'InternalModes method', modes.method);
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'InternalModes upper boundary', modes.upperBoundary);
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'InternalModes normalization', modes.normalization);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'number of particles or floats', nFloats);

% End definition mode
netcdf.endDef(ncid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Save output (internal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the data for the coordinate variables
netcdf.putVar(ncid, setprecision(iVarID), wavemodel.x);
netcdf.putVar(ncid, setprecision(jVarID), wavemodel.y);
netcdf.putVar(ncid, setprecision(kVarID), wavemodel.z); 

% Add the data for the dynamical variables
netcdf.putVar(ncid, setprecision(uVarID), [0 0 0 0], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
netcdf.putVar(ncid, setprecision(vVarID), [0 0 0 0], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
netcdf.putVar(ncid, setprecision(wVarID), [0 0 0 0], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
netcdf.putVar(ncid, setprecision(rhoprimeVarID), [0 0 0 0], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], rhoprime);
netcdf.putVar(ncid, setprecision(rhobarVarID), 0, wavemodel.Nz, rhobar);
netcdf.putVar(ncid, setprecision(drhobardzVarID), 0, wavemodel.Nz, drhobar_dz);
netcdf.putVar(ncid, setprecision(drhobardzdzVarID), 0, wavemodel.Nz, drhobar_dzdz);
% netcdf.putVar(ncid, setprecision(zetaVarID), [0 0 0 0], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
netcdf.putVar(ncid, setprecision(tVarID), 0, 1, time);
netcdf.putVar(ncid, setprecision(omegaForceVarID), 0, length(omega_forcing), omega_forcing);
netcdf.putVar(ncid, setprecision(ampForceVarID), 0, length(omega_forcing), amp_forcing);
netcdf.putVar(ncid, setprecision(kForceVarID), 0, length(omega_forcing), k);
netcdf.putVar(ncid, setprecision(lForceVarID), 0, length(omega_forcing), l);
netcdf.putVar(ncid, setprecision(phaseForceVarID), 0, length(omega_forcing), phase);
netcdf.putVar(ncid, setprecision(uForceVarID), [0 0], [wavemodel.Nz length(omega_forcing)], u_forcing);
netcdf.putVar(ncid, setprecision(vForceVarID), [0 0], [wavemodel.Nz length(omega_forcing)], v_forcing);
netcdf.putVar(ncid, setprecision(wForceVarID), [0 0], [wavemodel.Nz length(omega_forcing)], w_forcing);
netcdf.putVar(ncid, setprecision(rhoForceVarID), [0 0], [wavemodel.Nz length(omega_forcing)], rho_forcing);
netcdf.putVar(ncid, setprecision(xfloatVarID), 0, nFloats, x_float);
netcdf.putVar(ncid, setprecision(yfloatVarID), 0, nFloats, y_float);
netcdf.putVar(ncid, setprecision(zfloatVarID), 0, nFloats, z_float);

netcdf.close(ncid);	



