profilename = 'latmix-site1';
latitude = 31;
nOmega = 128;
nK = 128;

% m = 4;
nEVPMin = 2*256; % assumed minimum, can be overriden by the user
nEVPMax = 4*512;
nModes = 4*64;

[rho, N2, zIn] = InternalModes.StratificationProfileWithName(profilename);
z = linspace(min(zIn),max(zIn),5000);

im = InternalModesAdaptiveSpectral(rho,zIn,z,latitude,'nEVP',nEVPMax);
zInternal = im.z_xLobatto;
N2internal = im.N2_xLobatto;
N_max = max(sqrt(im.N2_xLobatto));
g = 9.81;

% Now re-inialize with the grid we extracted
nEVP = nEVPMin;
im = InternalModesAdaptiveSpectral(rho,zIn,zInternal,latitude,'nEVP',nEVP);
im.normalization = Normalization.kConstant;

% 
omega = exp(linspace(log(im.f0),log(0.99*N_max),nOmega));
k = zeros(1,nK);
k(2:nK) = exp(linspace(log(2*pi/1e7),log(1e1),nK-1));

methodName = 'ModesAtFrequency';
x = omega;
nX = length(x);

Phi = nan(length(zInternal),nX,nModes);
Gamma = nan(length(zInternal),nX,nModes);
h_out = nan(nX,nModes);

for i = 1:length(x)
    [F, G, h] = im.(methodName)(x(i));
    fprintf('%d of %d\n',i,length(x));
    % Increase the number of grid points until we get the
    % desired number of good quality modes (or reach some max).
    while( (isempty(h) || length(h) < nModes) && nEVP < nEVPMax )
        nEVP = nEVP + 128;
        im = InternalModesAdaptiveSpectral(rho,zIn,zInternal,latitude,'nEVP',nEVP);
        im.normalization = Normalization.kConstant;
        [F, G, h] = im.(methodName)(x(i));
    end
    if length(h) < nModes
        fprintf('Only found %d good modes (of %d requested). Proceeding anyway.\n',length(h),nModes);
    end
    j0 = min(length(h),nModes);
    
    h = reshape(h,1,[]);
    
    Phi(:,i,1:j0) = F(:,1:j0);
    Gamma(:,i,1:j0) = G(:,1:j0);
    h_out(i,1:j0)=h(1:j0);
end

F_omega = Phi;
G_omega = Gamma;
h_omega = h_out;
k_omega = sqrt(((omega.*omega - im.f0*im.f0).')./(g*h_omega));

methodName = 'ModesAtWavenumber';
x = k;
nX = length(x);

Phi = nan(length(zInternal),nX,nModes);
Gamma = nan(length(zInternal),nX,nModes);
h_out = nan(nX,nModes);

nEVP = nEVPMin;
for i = 1:length(x)
    [F, G, h] = im.(methodName)(x(i));
    fprintf('%d of %d\n',i,length(x));
    % Increase the number of grid points until we get the
    % desired number of good quality modes (or reach some max).
    while( (isempty(h) || length(h) < nModes) && nEVP < nEVPMax )
        nEVP = nEVP + 128;
        im = InternalModesAdaptiveSpectral(rho,zIn,zInternal,latitude,'nEVP',nEVP);
        im.normalization = Normalization.kConstant;
        [F, G, h] = im.(methodName)(x(i));
    end
    if length(h) < nModes
        fprintf('Only found %d good modes (of %d requested). Proceeding anyway.\n',length(h),nModes);
    end
    j0 = min(length(h),nModes);
    
    h = reshape(h,1,[]);
    
    Phi(:,i,1:j0) = F(:,1:j0);
    Gamma(:,i,1:j0) = G(:,1:j0);
    h_out(i,1:j0)=h(1:j0);
end

F_k = Phi;
G_k = Gamma;
h_k = h_out;
k2 = reshape(k .^2,[],1);
omega_k = sqrt(g * h_k .* k2 + im.f0*im.f0);

filename = sprintf('%s.mat',profilename);
save(filename,'F_omega','G_omega','omega', 'h_omega', 'k_omega', 'F_k','G_k','h_k', 'k', 'omega_k', 'latitude','zIn','N_max','zInternal','N2internal', 'rho', 'N2');