ReadOverNetwork = 0;
Nonlinear = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/';
end

if Nonlinear == 1
    runtype = 'NL';
    NonlinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart');
    model_dir = NonlinearSteadyStateFile;
else
    runtype = 'LIN';
    LinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_5xGM');
    model_dir = LinearSteadyStateFile;
end


% 
% WM = WintersModel(file);
% wavemodel = WM.wavemodel;

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
save(outputfile,'x','y','z','t','floatsPerLevel', 'model_dir');
