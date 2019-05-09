wintersOutputDirectory = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindown';

WM = WintersModel(wintersOutputDirectory);
wavemodel = WM.wavemodel;

nFiles = WM.NumberOf3DOutputFiles;
fileIncrements = 1:10:nFiles;

VortexEnergyConversionFactor = wavemodel.B_HKE_factor + wavemodel.B_PE_factor;
Vortex0EnergyConversionFactor = wavemodel.B0_HKE_factor;

nT = length(fileIncrements);
t = zeros(nT,1);
WaveEnergy_plus = zeros(nT,1);
WaveEnergy_minus = zeros(nT,1);
VortexEnergy = zeros(nT,1);
Vortex0Energy = zeros(nT,1);

startTime = datetime('now');
for iTime = 1:length(fileIncrements)
    if iTime>=2
        timePerStep = (datetime('now')-startTime)/(iTime-1);
        timeRemaining = (length(fileIncrements)-iTime+1)*timePerStep;
        fprintf('\treading values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(fileIncrements), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'dd:HH:MM:SS')) ;
    end
    if iTime == 2
        % The first time step takes extra long, because we're using a fixed
        % dimension length for time. So, let's reset the clock for
        % subsequent estimatation.
        startTime = datetime('now');
    end
    [ti,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(fileIncrements(iTime),'t','u','v','w','rho_prime');
    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(ti,u,v,rho_prime);
    
    Ap2 = wavemodel.Amp_plus .* conj(wavemodel.Amp_plus);
    Am2 = wavemodel.Amp_minus .* conj(wavemodel.Amp_minus);
    B2 = VortexEnergyConversionFactor .* (wavemodel.B .* conj(wavemodel.B));
    B02 = Vortex0EnergyConversionFactor.* (wavemodel.B0 .* conj(wavemodel.B0));
    
    t(iTime)=ti;
    WaveEnergy_plus(iTime) = sum(Ap2(:));
    WaveEnergy_minus(iTime) = sum(Am2(:));
    VortexEnergy(iTime) = sum(B2(:));
    Vortex0Energy(iTime) = sum(B02(:));
end

GM = GarrettMunkSpectrumConstantStratification(5.2e-3,[-4000 0], 33);
E_GM = GM.E;

figure
subplot(2,1,1)
plot(t/86400,WaveEnergy_plus), hold on
plot(t/86400,WaveEnergy_minus)
plot(t/86400,(E_GM/2)*ones(size(t)),'k')
xlim([min(t) max(t)]/86400)
ylim([0 1.05*max(max(WaveEnergy_minus),max(WaveEnergy_plus))])
legend('+\omega', '-\omega','E_{gm}/2', 'Location', 'Southeast')
ylabel('energy (m^3/s^2)')
title('Total depth integrated energy over time')

subplot(2,1,2)
plot(t/86400,VortexEnergy), hold on
plot(t/86400,Vortex0Energy)
xlim([min(t) max(t)]/86400)
legend('geostrophic', 'mean geostrophic')
xlabel('time (days)')
ylabel('energy (m^3/s^2)')

packfig(2,1)