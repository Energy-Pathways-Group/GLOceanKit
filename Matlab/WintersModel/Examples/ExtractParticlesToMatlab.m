baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';

% Version 2 files, from October 2018
LinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped');
NonlinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_NL_forced_damped');

model_dir = LinearSteadyStateFile;
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

outputfile = '/Users/jearly/Documents/ManuscriptRepositories/garrett-munk-lateral-diffusivity/data/2018_10/particles_LIN.mat';
save(outputfile,'x','y','z','t','floatsPerLevel', 'model_dir');
