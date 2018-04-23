% A few good plots:
% - show the decorrelation time as a function of wavenumber and mode.
% - show the wave energy fraction as a function of wavenumber and mode.
% - combine these somehow. AC is fraction of variance same as I.C. So, wave
% energy fraction, times AC. Then assess when it drops below, say 50%. This
% would tell you "How long linear IW's explain 50% of KE variance".



file = '/Volumes/seattle_data1/jearly/nsf_iwv/EarlyEtal_GM_NL_35e-11_36000s_restart_decomp.nc';
matfile = '/Volumes/seattle_data1/jearly/nsf_iwv/EarlyEtal_GM_NL_35e-11_36000s_restart_decomp.mat';

k = ncread(file, 'k');
l = ncread(file, 'l');
j = ncread(file, 'j');

[K,L,J] = ndgrid(k,l,j);
Kh = sqrt(K.*K + L.*L);

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(wavemodel.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

nK = length(kAxis)-1;
nModes = length(j);

t = ncread(file, 't');

variables = {'Ap_realp', 'Ap_imagp', 'Am_realp', 'Am_imagp', 'B_realp', 'B_imagp'};
Ppm_HKE_factor = wavemodel.Ppm_HKE_factor;
P0_HKE_factor = wavemodel.P0_HKE_factor;
conversion_factor = {Ppm_HKE_factor,Ppm_HKE_factor,Ppm_HKE_factor,Ppm_HKE_factor,P0_HKE_factor,P0_HKE_factor};
Nvars = length(variables);

% Not sure what a good value to use is.
decorrelationCutoff = 0.5;

waveHKE = zeros(nK,nModes);
vortexHKE = zeros(nK,nModes);
waveDecorrelationTime = zeros(nK,nModes);
vortexDecorrelationTime = zeros(nK,nModes);
% for iK = 1:nK
for iMode = 1:1
    fprintf('iMode: %d\n',iMode);
    for iK = 1:1:nK
        fprintf('\tiK: %d\n',iK);
        indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1) );
        
        AC = zeros(length(t),Nvars);
        DOF = zeros(length(t),Nvars);
        nloops = zeros(1,Nvars);
        HKE = zeros(1,Nvars);
        for iIndex = 1:length(indicesForK)
            [i,j] = ind2sub([size(K,1) size(K,2)],indicesForK(iIndex));
            
            for iVar = 1:Nvars
                u = double(squeeze(ncread(file, variables{iVar}, [i j iMode 1], [1 1 1 length(t)], [1 1 1 1])));
                [ACu, DOFu] = Autocorrelation(u, length(t)-1);
                if any(isnan(ACu))
                    continue; % this will occur for the occasional unresolved mode. Seems to only be the Nyquist, which is okay.
                end
                AC(:,iVar) = AC(:,iVar) + ACu;
                DOF(:,iVar) = DOF(:,iVar) + DOFu;
                nloops(1,iVar) = nloops(1,iVar)+1;
                
                HKE(iVar) = HKE(iVar) + mean(u.*conj(u))*conversion_factor{iVar}(i,j,iMode);
            end
            
        end
        AC = AC./nloops;
        waveAC = mean(AC(:,1:4),2);
        waveDOF = sum(DOF(:,1:4),2);
        waveSE =  circshift(sqrt((1 + 2*cumsum(waveAC.^2,1))./waveDOF),1);
        waveSE(2) = sqrt(1./DOF(1));
        waveSE(1) = 0;
        
        i = find( waveAC < decorrelationCutoff & waveAC > waveSE,1,'first');
        if isempty(i)
            waveDecorrelationTime(iK,iMode) = Inf;
        else
            waveDecorrelationTime(iK,iMode) = t(i);
        end
        
        vortexAC = mean(AC(:,5:6),2);
        vortexDOF = sum(DOF(:,5:6),2);
        vortexSE =  circshift(sqrt((1 + 2*cumsum(vortexAC.^2,1))./vortexDOF),1);
        vortexSE(2) = sqrt(1./DOF(1));
        vortexSE(1) = 0;
        
        i = find( vortexAC < decorrelationCutoff & vortexAC > vortexSE,1,'first');
        if isempty(i)
            vortexDecorrelationTime(iK,iMode) = Inf;
        else
            vortexDecorrelationTime(iK,iMode) = t(i);
        end
        
        waveHKE(iK,iMode) = sum(HKE(1:4));
        vortexHKE(iK,iMode) = sum(HKE(5:6));
    end
end

save(matfile,'waveHKE','vortexHKE','waveDecorrelationTime','vortexDecorrelationTime');

HKE_fraction = waveHKE./(waveHKE+vortexHKE);

lengthScaleAxis = 2*pi./kAxis(1:nK)/1000;

figure
subplot(2,1,1)
plot(lengthScaleAxis,squeeze(HKE_fraction(:,iMode)), 'LineWidth',1,'Color',0.0*[1.0 1.0 1.0]), xlog
xlabel('km')
title('wave energy fraction')

subplot(2,1,2)
timescale = 1/86400;
plot(lengthScaleAxis,timescale*squeeze(waveDecorrelationTime(:,iMode)), 'LineWidth',1.5), hold on
plot(lengthScaleAxis,timescale*squeeze(vortexDecorrelationTime(:,iMode)), 'LineWidth',1.5), xlog
xlabel('km')
ylabel('days')
title('linear solution decorrelation time')
legend('wave','vortex')
