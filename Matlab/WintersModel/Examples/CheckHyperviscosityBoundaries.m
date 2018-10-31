scaleFactor = 1;
LoadFigureDefaults;

runtype = 'linear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/model_raw/';
end

if strcmp(runtype,'linear')
    file = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped');
elseif strcmp(runtype,'nonlinear')
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped'); 
else
    error('invalid run type.');
end

WM = WintersModel(file);
wavemodel = WM.wavemodel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Do the wave-vortex decomposition for t=0
%

% These first two lines actually do the decomposition
[t0,u,v,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t0,u,v,rho_prime);

waveHKE_full = wavemodel.Ppm_HKE_factor .* ( abs(wavemodel.Amp_plus).^2 + abs(wavemodel.Amp_minus).^2 );
waveVKE_full = wavemodel.Ppm_VKE_factor .* ( abs(wavemodel.Amp_plus).^2 + abs(wavemodel.Amp_minus).^2 );
vortexHKE_full = wavemodel.P0_HKE_factor .* abs(wavemodel.B).^2;

% Now we need to figure out how to display this information...

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(wavemodel.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

% Thi is the final output axis for wavenumber
k = reshape(kAxis(1:(length(kAxis)-1)),[],1);

% Mode axis is just what we already have
j = wavemodel.j;

Kh = wavemodel.Kh;
RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);

nK = length(k);
nModes = length(j);

waveHKE_t0 = zeros(nK,nModes);
waveVKE_t0 = zeros(nK,nModes);
vortexHKE_t0 = zeros(nK,nModes);
for iMode = 1:1:nModes
    for iK = 1:1:nK
        indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1)  & ~squeeze(RedundantCoefficients(:,:,1)) );
        for iIndex = 1:length(indicesForK)
            [i,m] = ind2sub([size(Kh,1) size(Kh,2)],indicesForK(iIndex));
            waveHKE_t0(iK,iMode) = waveHKE_t0(iK,iMode) + waveHKE_full(i,m,iMode);
            waveVKE_t0(iK,iMode) = waveVKE_t0(iK,iMode) + waveVKE_full(i,m,iMode);
            vortexHKE_t0(iK,iMode) = vortexHKE_t0(iK,iMode) + vortexHKE_full(i,m,iMode);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Do the wave-vortex decomposition for t=t_max
%

maxFileIndex = WM.NumberOf3DOutputFiles;

% These first two lines actually do the decomposition
[t_max,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(maxFileIndex,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t_max,u,v,rho_prime);

waveHKE_full = wavemodel.Ppm_HKE_factor .* ( abs(wavemodel.Amp_plus).^2 + abs(wavemodel.Amp_minus).^2 );
waveVKE_full = wavemodel.Ppm_VKE_factor .* ( abs(wavemodel.Amp_plus).^2 + abs(wavemodel.Amp_minus).^2 );
vortexHKE_full = wavemodel.P0_HKE_factor .* abs(wavemodel.B).^2;

waveHKE = zeros(nK,nModes);
waveVKE = zeros(nK,nModes);
vortexHKE = zeros(nK,nModes);
for iMode = 1:1:nModes
    for iK = 1:1:nK
        indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1)  & ~squeeze(RedundantCoefficients(:,:,1)) );
        for iIndex = 1:length(indicesForK)
            [i,m] = ind2sub([size(Kh,1) size(Kh,2)],indicesForK(iIndex));
            waveHKE(iK,iMode) = waveHKE(iK,iMode) + waveHKE_full(i,m,iMode);
            waveVKE(iK,iMode) = waveVKE(iK,iMode) + waveVKE_full(i,m,iMode);
            vortexHKE(iK,iMode) = vortexHKE(iK,iMode) + vortexHKE_full(i,m,iMode);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the damping scales
%

% Kraig Winters uses an e-fold time to set \nu_x in the hypervicous
% operator. We start by computing \nu_x and \nu_z.
dx = wavemodel.x(2)-wavemodel.x(1);
dz = wavemodel.z(2)-wavemodel.z(1);
nu_x = (-1)^(WM.p_x+1)*power(dx/pi,2*WM.p_x) / WM.T_diss_x;
nu_z = (-1)^(WM.p_z+1)*power(dz/pi,2*WM.p_z) / WM.T_diss_z;

threshold = 0.6; % exp(-1) should get you all the way to the Nyquist.
tau = max(WM.T_diss_x,WM.T_diss_z);
m = (pi/(length(wavemodel.j)*dz))*wavemodel.j;
lambda_x = nu_x*(sqrt(-1)*wavemodel.k).^(2*WM.p_x);
lambda_z = nu_z*(sqrt(-1)*m).^(2*WM.p_z);
maxKindex = find(exp(lambda_x*tau)<threshold,1,'first');
maxJindex = find(exp(lambda_z*tau)<threshold,1,'first');

maxK = wavemodel.k(maxKindex);
maxMode = maxJindex;

nK = length(wavemodel.k)/2 + 1;
k_diss = abs(wavemodel.k(1:nK));
j_diss = 0:max(wavemodel.j); % start at 0, to help with contour drawing
[K,J] = ndgrid(k_diss,j_diss);
M = (2*pi/(length(wavemodel.j)*dz))*J/2;

lambda_x = nu_x*(sqrt(-1)*K).^(2*WM.p_x);
lambda_z = nu_z*(sqrt(-1)*M).^(2*WM.p_z);
% tau = WM.VariableFieldsFrom3DOutputFileAtIndex(WM.NumberOf3DOutputFiles,'t');
tau = max(WM.T_diss_x,WM.T_diss_z);
R = exp(2*(lambda_x+lambda_z)*(t_max-t0));

% The highest wavenumber should e-fold in time tau, so let's contour the
% area that retains 90% of its value
C = contourc(j_diss,k_diss,R,0.1*[1 1]);
n = C(2,1);
j_damp = C(1,1+1:n);
k_damp = C(2,1+1:n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the axis labels
%

ticks_x = [100;10];
labels_x = cell(length(ticks_x),1);
for i=1:length(ticks_x)
    labels_x{i} = sprintf('%d',round(ticks_x(i)));
end
ticks_x = 2*pi./(1e3*ticks_x);

ticks_y = [1;10;100];
labels_y = cell(length(ticks_y),1);
for i=1:length(ticks_y)
    labels_y{i} = sprintf('%d',round(ticks_y(i)));
end

figure
subplot(1,3,1)
pcolor( k, j, log10((waveHKE./waveHKE_t0).') ), xlog, ylog, shading flat, hold on
caxis([-3 0])
plot( k_damp, j_damp, 'LineWidth', 4, 'Color', 0*[1 1 1])
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

subplot(1,3,2)
pcolor( k, j, log10((waveVKE./waveVKE_t0).') ), xlog, ylog, shading flat, hold on
caxis([-3 0])
plot( k_damp, j_damp, 'LineWidth', 4, 'Color', 0*[1 1 1])
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

subplot(1,3,3)
pcolor( k, j, log10((waveHKE./waveHKE_t0).')+log10((waveVKE./waveVKE_t0).') ), xlog, ylog, shading flat, hold on
caxis([-3 0])
plot( k_damp, j_damp, 'LineWidth', 4, 'Color', 0*[1 1 1])
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

