%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SaveModelOutputToNetCDF
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% December 6th, 2016      Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 800e3;
Ly = 100e3;
Lz = 5000;

Nx = 2048/2;
Ny = 256/2;
Nz = 128;

% Lx = 15e3;
% Ly = 15e3;
% Lz = 5000;
% 
% Nx = 32;
% Ny = 32;
% Nz = 16;

latitude = 31;
N0 = 5.2e-3;
GMReferenceLevel = 1.0;

timeStep = 15*60; % in seconds
maxTime = 4*86400;

outputfolder = '/Volumes/OceanTransfer';
outputfolder = '/Users/jearly/Desktop';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.ShowDiagnostics();
wavemodel.InitializeWithGMSpectrum(GMReferenceLevel);

t = (0:timeStep:maxTime)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/InternalWaveModel_%dx%dx%d@%s.nc', outputfolder,Nx,Ny,Nz,datestr(datetime('now'),'yyyy-mm-dd_HH:MM:SS'));

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 4;
totalSize = totalFields*bytePerFloat*length(t)*(wavemodel.Nx)*(wavemodel.Ny)*(wavemodel.Nz)/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',filepath,totalSize);

ncid = netcdf.create(filepath, 'CLOBBER');

% Define the dimensions
xDimID = netcdf.defDim(ncid, 'x', wavemodel.Nx);
yDimID = netcdf.defDim(ncid, 'y', wavemodel.Ny);
zDimID = netcdf.defDim(ncid, 'z', wavemodel.Nz);
tDimID = netcdf.defDim(ncid, 't', netcdf.getConstant('NC_UNLIMITED'));

% Define the coordinate variables
xVarID = netcdf.defVar(ncid, 'x', ncPrecision, xDimID);
yVarID = netcdf.defVar(ncid, 'y', ncPrecision, yDimID);
zVarID = netcdf.defVar(ncid, 'z', ncPrecision, zDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
netcdf.putAtt(ncid,xVarID, 'units', 'm');
netcdf.putAtt(ncid,yVarID, 'units', 'm');
netcdf.putAtt(ncid,zVarID, 'units', 'm');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the dynamical variables
uVarID = netcdf.defVar(ncid, 'u', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
vVarID = netcdf.defVar(ncid, 'v', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
wVarID = netcdf.defVar(ncid, 'w', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
zetaVarID = netcdf.defVar(ncid, 'zeta', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
netcdf.putAtt(ncid,uVarID, 'units', 'm/s');
netcdf.putAtt(ncid,vVarID, 'units', 'm/s');
netcdf.putAtt(ncid,wVarID, 'units', 'm/s');
netcdf.putAtt(ncid,zetaVarID, 'units', 'm');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', N0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

% End definition mode
netcdf.endDef(ncid);

% Add the data for the coordinate variables
netcdf.putVar(ncid, setprecision(xVarID), wavemodel.x);
netcdf.putVar(ncid, setprecision(yVarID), wavemodel.y);
netcdf.putVar(ncid, setprecision(zVarID), wavemodel.z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the model, and write the output to NetCDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = datetime('now');
fprintf('Starting numerical simulation on %s\n', datestr(startTime));
for iTime=1:length(t)
    if mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-1);
        timeRemaining = (length(t)-iTime+1)*timePerStep;   
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end
    [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
    [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
    netcdf.putVar(ncid, setprecision(uVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
    netcdf.putVar(ncid, setprecision(vVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
    netcdf.putVar(ncid, setprecision(wVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
    netcdf.putVar(ncid, setprecision(zetaVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
    netcdf.putVar(ncid, setprecision(tVarID), iTime-1, 1, t(iTime));
end
fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
