% A few good plots:
% - show the decorrelation time as a function of wavenumber and mode.
% - show the wave energy fraction as a function of wavenumber and mode.
% - combine these somehow. AC is fraction of variance same as I.C. So, wave
% energy fraction, times AC. Then assess when it drops below, say 50%. This
% would tell you "How long linear IW's explain 50% of KE variance".



file = '/Volumes/seattle_data1/jearly/nsf_iwv/EarlyEtal_GM_NL_35e-11_36000s_restart_decomp.nc';

k = ncread(file, 'k');
l = ncread(file, 'l');
j = ncread(file, 'j');

[K,L,J] = ndgrid(k,l,j);
Kh = sqrt(K.*K + L.*L);

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(wavemodel.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

t = ncread(file, 't');

variables = {'Ap_realp', 'Ap_imagp', 'Am_realp', 'Am_imagp', 'B_realp', 'B_imagp'};
Ppm_HKE_factor = wavemodel.Ppm_HKE_factor;
P0_HKE_factor = wavemodel.P0_HKE_factor;
conversion_factor = {Ppm_HKE_factor,Ppm_HKE_factor,Ppm_HKE_factor,Ppm_HKE_factor,P0_HKE_factor,P0_HKE_factor};
Nvars = length(variables);

AC = zeros(length(t),Nvars);
DOF = zeros(length(t),Nvars);
HKE = zeros(Nvars,1);
nloops = zeros(1,Nvars);
% for iK = 1:length(kAxis)
for iMode = 1:1
    for iK = 50:50
        indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1) );
        
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
    end
end

AC = AC./nloops;

% let's combine the real and imaginary parts
for i=Nvars:-2:1
    AC(:,i-1) = (AC(:,i) + AC(:,i-1))/2;
    DOF(:,i-1) = DOF(:,i) + DOF(:,i-1);
    AC(:,i) = [];
    DOF(:,i) = [];
    HKE(i-1) = HKE(i) + HKE(i-1);
    HKE(i) = [];
end

SE_indep = t(2:end);
SE =  sqrt((1 + 2*cumsum(AC.^2,1))./DOF);
SE(1,:) = sqrt(1./DOF(1,:)); % first point is lag 1
SE(end,:) = []; % there is no end point

figure
plot(t,AC, 'LineWidth',1)
hold on
plot(SE_indep, [3*SE,-3*SE], 'LineWidth', 1.5, 'Color',0.4*[1.0 1.0 1.0] )
xlabel('time lag (seconds)')
ylabel('temporal correlation')
legend('A_p', 'A_m', 'B', 'SE_p', 'SE_m', 'SE_0')
title(sprintf('%d km wavelength, %d%% wave energy', round(2*pi/kAxis(iK)/1000), round(100*(HKE(1) + HKE(2))/sum(HKE))) )