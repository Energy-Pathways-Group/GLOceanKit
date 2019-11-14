% % Fetch a pre-built profile from the internal modes class
% [rho, N2, zIn] = InternalModes.StratificationProfileWithName('exponential');
% latitude = 33;
% 
% GM = GarrettMunkSpectrum('/Users/jearly/Documents/ProjectRepositories/GLOceanKit/Matlab/InternalWaveSpectrum/PrecomputedProfiles/exponential-free-surface');
% 
% % Initialize the GarrettMunkSpectrum class with the profile, but do *not*
% % reinitialize it if it already exists. Recomputing the modes can take a
% % while.
% if ~exist('GM','var')
%     GM = GarrettMunkSpectrum(rho,zIn,latitude);
% end

omega = linspace(0,GM.N_max,2000).';
[S,k] = GM.SSHSpectrumAtFrequencies(omega);
figure, plot(omega,sum(S,2,'omitnan')), ylog, xlog

S_zeta = GM.IsopycnalSpectrumAtFrequencies(0,omega);
hold on, plot(omega,S_zeta)

k_edges = exp(linspace(log(2*pi/1e10),log(2*pi/15e1),101));
omega_edges = exp(linspace(log(GM.f0/10),log(GM.N_max),101));
k_axis = k_edges(1:end-1)+diff(k_edges);
omega_axis = omega_edges(1:end-1)+diff(omega_edges);

S_omega_k = zeros(length(omega_axis),length(k_axis));
omega_omega_k = repmat(omega,1,size(k,2));

for i=1:length(omega_axis)
    for j=1:length(k_axis)
        omega_indices = omega_omega_k > omega_edges(i) & omega_omega_k <= omega_edges(i+1);
        k_indices = k > k_edges(j) & k <= k_edges(j+1);
        S_omega_k(i,j) = sum(S(omega_indices & k_indices))/(k_edges(j+1) - k_edges(j));
    end
end


k_flat = k(:);
[N,edges,bin] = histcounts(k_flat,k_edges);