file = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-3.nc';

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

figure
for i=1:10
    netcdfTools.SetWaveModelToIndex(i);
    [k,j,A0KJ] = wvm.ConvertToWavenumberAndMode(abs(wvm.A0).^2);
    plot(k,sum(A0KJ,2)), ylog, hold on
end