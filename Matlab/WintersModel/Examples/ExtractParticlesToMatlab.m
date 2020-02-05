ReadOverNetwork = 1;
Nonlinear = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/';
end

if Nonlinear == 1
    runtype = 'NL';
    NonlinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_03xGM');
    model_dir = NonlinearSteadyStateFile;
else
    runtype = 'LIN';
    LinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_fresh');
    model_dir = LinearSteadyStateFile;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Comput the actual energy of the internal wave field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WM = WintersModel(file);
wavemodel = WM.wavemodel;
[t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
GMEnergyLevel = sum( abs(wavemodel.Amp_minus(:)).^2 + abs(wavemodel.Amp_plus(:)).^2 )/E;
fprintf('Measured energy level: %.2f GM\n', GMEnergyLevel);


eulerian_file = [model_dir 'input/SaveIC_EarlyIWmodel.nc'];
lagrangian_dir = [model_dir '/output/lagrangian/'];
floatsPerLevel = 100;

UniqueParticleFiles = dir([lagrangian_dir 'particles_*.nc']);

for iFile = 1:length(UniqueParticleFiles)
    file = [lagrangian_dir UniqueParticleFiles(iFile).name];
    if (iFile == 1)
        t = ncread(file,'t_secs');
        x = ncread(file,'x')';
        y = ncread(file,'y')';
        z = ncread(file,'z')';
    else
        x = cat(2,x,ncread(file,'x')');
        y = cat(2,y,ncread(file,'y')');
        z = cat(2,z,ncread(file,'z')');
    end
end

outputfile = sprintf('%s_particles.mat',model_dir);
save(outputfile,'x','y','z','t','floatsPerLevel', 'model_dir','GMEnergyLevel');
