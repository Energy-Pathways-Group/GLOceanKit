ReadOverNetwork = 1;
Nonlinear = 1;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2019_12/';
end

if Nonlinear == 1
    runtype = 'NL';
    NonlinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_01xGM');
    model_dir = NonlinearSteadyStateFile;
else
    runtype = 'LIN';
    LinearSteadyStateFile = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_01xGM');
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

outputfile = sprintf('/Users/jearly/Documents/ManuscriptRepositories/garrett-munk-lateral-diffusivity/data/2020_01/particles_%s.mat',runtype);
save(outputfile,'x','y','z','t','floatsPerLevel', 'model_dir');
