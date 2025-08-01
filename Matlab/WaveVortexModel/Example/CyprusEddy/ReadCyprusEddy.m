file = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-more-stratification.nc';

netcdfTools = WaveVortexModelNetCDFTools(file);
wvm = netcdfTools.InitializeWaveVortexModelFromNetCDFFile();

% netcdfTools.SetWaveModelToIndex(3);

% [k,j,ApKJ,AmKJ,A0KJ] = wvm.ConvertToWavenumberAndMode(abs(wvm.Ap).^2,abs(wvm.Am).^2,abs(wvm.A0).^2);
% 
% figure
% subplot(2,2,1)
% plot(k,sum(ApKJ,2)), hold on, plot(k,sum(AmKJ,2)), ylog
% subplot(2,2,2)
% plot(k,sum(A0KJ,2)), ylog
% 
% subplot(2,2,3)
% plot(j,sum(ApKJ,1)), hold on, plot(j,sum(AmKJ,1)), ylog
% subplot(2,2,4)
% plot(j,sum(A0KJ,1)), ylog

% [Ep,Em,E0] = wvm.energyFlux;
% inertialFlux = sum(Ep(1,1,:)) + sum(Em(1,1,:));
% Ep(1,1,:) = 0; Em(1,1,:) = 0;
% waveFlux = sum(Ep(:)) + sum(Em(:));
% fprintf('total spectral: %g, total depth integrated: %g\n',wvm.totalEnergy,wvm.totalHydrostaticEnergy);
% fprintf('io flux: %g, wave flux: %g, geostrophic flux: %g. Net: %g\n',inertialFlux,waveFlux,sum(E0(:)),inertialFlux+waveFlux+sum(E0(:)))

figure
sp1 = subplot(1,2,1);
sp2 = subplot(1,2,2);
for i=1:10
    netcdfTools.SetWaveModelToIndex(i);
  % [k,j,A0KJ] = wvm.ConvertToWavenumberAndMode(abs(wvm.Ap).^2+abs(wvm.Am).^2);
  %[k,j,A0KJ] = wvm.ConvertToWavenumberAndMode(wvm.A0_TE_factor.*abs(wvm.A0).^2);
  [Ep,Em,E0] = wvm.energyFlux;
    [k,j,A0KJ] = wvm.ConvertToWavenumberAndMode(E0);
    subplot(sp1);
    plot(j,sum(A0KJ,1).'),xlog, hold on
    subplot(sp2);
    plot(k,sum(A0KJ,2)),xlog, hold on
     pause(1)
end

% xlabel('mode number')
xlabel('horizontal wavenumber')
ylabel('geostrophic energy')

return

[k,j,ApKJ,AmKJ,A0KJ] = wvm.ConvertToWavenumberAndMode(abs(wvm.Ap).^2,abs(wvm.Am).^2,abs(wvm.A0).^2);
% [k,j,ApKJ,AmKJ,A0KJ] = self.ConvertToWavenumberAndMode(abs(uNLbar.^2),abs(vNLbar).^2,abs(nNLbar).^2);
[k,j,ApKJ,AmKJ,A0KJ] = wvm.ConvertToWavenumberAndMode(Ep,Em,E0);
figure
subplot(2,2,1)
plot(k,sum(ApKJ,2)), hold on, plot(k,sum(AmKJ,2))%, ylog
subplot(2,2,2)
plot(k,sum(A0KJ,2))%, ylog

subplot(2,2,3)
plot(j,sum(ApKJ,1)), hold on, plot(j,sum(AmKJ,1))%, ylog
subplot(2,2,4)
plot(j,sum(A0KJ,1))%, ylog