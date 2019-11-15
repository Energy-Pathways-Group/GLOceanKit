% Fetch a pre-built profile from the internal modes class
[rho, N2, zIn] = InternalModes.StratificationProfileWithName('exponential');
latitude = 33;

GM = GarrettMunkSpectrum('/Users/jearly/Documents/ProjectRepositories/GLOceanKit/Matlab/InternalWaveSpectrum/PrecomputedProfiles/exponential-free-surface');

% Initialize the GarrettMunkSpectrum class with the profile, but do *not*
% reinitialize it if it already exists. Recomputing the modes can take a
% while.
if ~exist('GM','var')
    GM = GarrettMunkSpectrum(rho,zIn,latitude);
end

% omega = linspace(0,GM.N_max,2000).';
% [S,k] = GM.SSHSpectrumAtFrequencies(omega);
% figure, plot(omega,sum(S,2,'omitnan')), ylog, xlog
% 
% S_zeta = GM.IsopycnalSpectrumAtFrequencies(0,omega);
% hold on, plot(omega,S_zeta)
% 
% k_edges = exp(linspace(log(2*pi/1e10),log(2*pi/15e1),101));
% omega_edges = exp(linspace(log(GM.f0/10),log(GM.N_max),101));
% k_axis = k_edges(1:end-1)+diff(k_edges);
% omega_axis = omega_edges(1:end-1)+diff(omega_edges);
% 
% S_omega_k = zeros(length(omega_axis),length(k_axis));
% omega_omega_k = repmat(omega,1,size(k,2));
% 
% for i=1:length(omega_axis)
%     for j=1:length(k_axis)
%         omega_indices = omega_omega_k > omega_edges(i) & omega_omega_k <= omega_edges(i+1);
%         k_indices = k > k_edges(j) & k <= k_edges(j+1);
%         S_omega_k(i,j) = sum(S(omega_indices & k_indices))/(k_edges(j+1) - k_edges(j));
%     end
% end
% 
% 
% 
% figure, pcolor(k_axis,omega_axis,log(S_omega_k.')), shading flat,  xlog, ylog
% xlim([min(k_axis) max(k_axis)])
% ylim([min(omega_axis) max(omega_axis)])

[Sssh_kj, k, j] = GM.SSHSpectrumAtFrequencyWavenumbers;

figure, pcolor(k,j,log(Sssh_kj).'), shading flat, xlog

[Sssh_k_omega, k_, omega_] = GM.SSHSpectrumAtWavenumberAndFrequency();
trapz(k_,trapz(omega_,Sssh_k_omega,2),1)


% figure, pcolor(k_,omega_,log(Sssh_k_omega).'), shading flat, xlog, ylog
figure, pcolor(k_,omega_,(log(Sssh_k_omega.*(omega_).').') ), shading flat, xlog, ylog
xlim([min(1e3/1e7) max(1e3/15e3)])
xticks([1e-3 1e-2 1/16])
labels = cell(3,1); labels{1} = '1000 km';  labels{2} = '100 km';  labels{3} = '15 km';
xticklabels(labels)

ticks = exp(linspace(log(GM.f0),log(GM.N_max),5));
labels = MakeLabelsForFrequencies(ticks,GM.f0,GM.N_max);
yticks(ticks)
yticklabels(labels)

cmap = colormap();
cmap(1,:) = 1;
colormap(cmap)
