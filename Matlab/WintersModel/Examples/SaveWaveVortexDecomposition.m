NonlinearSpindownFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
NonlinearForcedFromInitialConditionsFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
LinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';
NonlinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';

file = LinearSteadyStateFile;

output_directory = '/Volumes/seattle_data1/jearly/nsf_iwv';
output_directory = '/Users/jearly/Documents';

shouldChunk = 0;
[filepath,name,ext] = fileparts(file);
outputfile = fullfile(output_directory,strcat(name,'_decomp.nc'));

% Notes on chunking, while reading temporal slices
%   Network read, chunked:         144.6s
%   Network read, not chunked:     37.9s
%   Local read (SSD), chunked:     5.6s
%   Local read (SSD), not chunked: 2.5s
% So, I'm not seeing an advantage with the sizes that I chose.

WM = WintersModel(file);
wavemodel = WM.WaveModelFromInitialConditionsFile;

nFiles = WM.NumberOfTimeSteps;
fileIncrements = 1:2:nFiles;
% fileIncrements = 1:2;

Nk = length(wavemodel.k);
Nl = length(wavemodel.l);
Nj = length(wavemodel.j);
Nt = length(fileIncrements);

precision = 'single';

if strcmp(precision,'single')
    ncPrecision = 'NC_FLOAT';
    setprecision = @(x) single(x);
    bytePerFloat = 4;
else
    ncPrecision = 'NC_DOUBLE';
    setprecision = @(x) double(x);
    bytePerFloat = 8;
end

% Chunking: https://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
[csize, nelems, premp] = netcdf.getChunkCache();
D = csize/bytePerFloat;
c = (D/(Nk*Nl*Nj*Nt))^(1/4);
chunkSize = floor(c*[Nk Nl Nj Nt]);

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(outputfile, cmode);

% Define the dimensions
kDimID = netcdf.defDim(ncid, 'k', Nk);
lDimID = netcdf.defDim(ncid, 'l', Nl);
jDimID = netcdf.defDim(ncid, 'j', Nj);
tDimID = netcdf.defDim(ncid, 't', Nt);

% Define the coordinate variables
kVarID = netcdf.defVar(ncid, 'k', ncPrecision, kDimID);
lVarID = netcdf.defVar(ncid, 'l', ncPrecision, lDimID);
jVarID = netcdf.defVar(ncid, 'j', ncPrecision, jDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
netcdf.putAtt(ncid,kVarID, 'units', 'radians/m');
netcdf.putAtt(ncid,lVarID, 'units', 'radians/m');
netcdf.putAtt(ncid,jVarID, 'units', 'mode number');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the wave-vortex variables
ApRealVarID = netcdf.defVar(ncid, 'Ap_realp', ncPrecision, [kDimID,lDimID,jDimID,tDimID]);
ApImagVarID = netcdf.defVar(ncid, 'Ap_imagp', ncPrecision, [kDimID,lDimID,jDimID,tDimID]);
AmRealVarID = netcdf.defVar(ncid, 'Am_realp', ncPrecision, [kDimID,lDimID,jDimID,tDimID]);
AmImagVarID = netcdf.defVar(ncid, 'Am_imagp', ncPrecision, [kDimID,lDimID,jDimID,tDimID]);
BRealVarID = netcdf.defVar(ncid, 'B_realp', ncPrecision, [kDimID,lDimID,jDimID,tDimID]);
BImagVarID = netcdf.defVar(ncid, 'B_imagp', ncPrecision, [kDimID,lDimID,jDimID,tDimID]);
if shouldChunk == 1
    netcdf.defVarChunking(ncid,ApRealVarID,'CHUNKED',chunkSize);
    netcdf.defVarChunking(ncid,ApImagVarID,'CHUNKED',chunkSize);
    netcdf.defVarChunking(ncid,AmRealVarID,'CHUNKED',chunkSize);
    netcdf.defVarChunking(ncid,AmImagVarID,'CHUNKED',chunkSize);
    netcdf.defVarChunking(ncid,BRealVarID,'CHUNKED',chunkSize);
    netcdf.defVarChunking(ncid,BImagVarID,'CHUNKED',chunkSize);
end
netcdf.putAtt(ncid,ApRealVarID, 'units', 'm/s');
netcdf.putAtt(ncid,ApImagVarID, 'units', 'm/s');
netcdf.putAtt(ncid,AmRealVarID, 'units', 'm/s');
netcdf.putAtt(ncid,AmImagVarID, 'units', 'm/s');
netcdf.putAtt(ncid,BRealVarID, 'units', 'm');
netcdf.putAtt(ncid,BImagVarID, 'units', 'm');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', wavemodel.latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', wavemodel.N0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Lz', wavemodel.Lz);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'source-file', file);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

% End definition mode
netcdf.endDef(ncid);

% Add the data for the coordinate variables
netcdf.putVar(ncid, kDimID, wavemodel.k);
netcdf.putVar(ncid, lDimID, wavemodel.l);
netcdf.putVar(ncid, jDimID, wavemodel.j);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 6;
totalSize = totalFields*bytePerFloat*length(fileIncrements)*Nk*Nl*Nj/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',outputfile,totalSize);

startTime = datetime('now');
for iTime = 1:length(fileIncrements)
    if iTime>=2 %|| mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-1);
        timeRemaining = (length(fileIncrements)-iTime+1)*timePerStep;
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(fileIncrements), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'dd:HH:MM:SS')) ;
    end
    if iTime == 2
        % The first time step takes extra long, because we're using a fixed
        % dimension length for time. So, let's reset the clock for
        % subsequent estimatation.
        startTime = datetime('now');
    end
    [t,u,v,w,rho_prime] = WM.VariableFieldsAtTimeIndex(fileIncrements(iTime),'t','u','v','w','rho_prime');

    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
    
    netcdf.putVar(ncid, ApRealVarID, [0 0 0 iTime-1], [Nk Nl Nj 1], real(wavemodel.Amp_plus));
    netcdf.putVar(ncid, ApImagVarID, [0 0 0 iTime-1], [Nk Nl Nj 1], imag(wavemodel.Amp_plus));
    netcdf.putVar(ncid, AmRealVarID, [0 0 0 iTime-1], [Nk Nl Nj 1], real(wavemodel.Amp_minus));
    netcdf.putVar(ncid, AmImagVarID, [0 0 0 iTime-1], [Nk Nl Nj 1], imag(wavemodel.Amp_minus));
    netcdf.putVar(ncid, BRealVarID, [0 0 0 iTime-1], [Nk Nl Nj 1], real(wavemodel.B));
    netcdf.putVar(ncid, BImagVarID, [0 0 0 iTime-1], [Nk Nl Nj 1], imag(wavemodel.B));
    
    netcdf.putVar(ncid, tVarID, iTime-1, 1, t);
end

netcdf.close(ncid);	