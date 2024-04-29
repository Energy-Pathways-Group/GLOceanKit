function wvt = waveVortexTransformFromFile(path,options)
% Initialize a WVTransformHydrostatic instance from an existing file
%
% This static method is called by WVTransform.waveVortexTransformFromFile
% and should not need to be called directly.
%
% - Topic: Initialization (Static)
% - Declaration: wvt = waveVortexTransformFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter iTime: (optional) time index to initialize from (default 1)
arguments
    path char {mustBeFile}
    options.iTime (1,1) double {mustBePositive} = 1
end

ncfile = NetCDFFile(path);

requiredVariables = {'x','y','z','Lx','Ly','Lz','j','latitude','t0','t','rho0','dLnN2','PFinv','QGinv','PF','QG','P','Q','h','shouldAntialias'};
if ~all(isKey(ncfile.variableWithName,requiredVariables) | isKey(ncfile.attributes,requiredVariables))
    error('This files is missing required variables or attributes to initialize a WVTransform.')
end

[x,y,z,Lx,Ly,Lz,rho0,latitude,dLnN2,PFinv,QGinv,PF,QG,P,Q,h,shouldAntialias] = ncfile.readVariables('x','y','z','Lx','Ly','Lz','rho0','latitude','dLnN2','PFinv','QGinv','PF','QG','P','Q','h','shouldAntialias');
Nx = length(x);
Ny = length(y);
Nz = length(z);

[filepath,name,~] = fileparts(path);
if isempty(filepath)
    matFilePath = sprintf('%s.mat',name);
else
    matFilePath = sprintf('%s/%s.mat',filepath,name);
end
if ~isfile(matFilePath)
    error('The .mat sidecar file is missing, which is necessary for the hydrostatic transformations.')
end
matFile = load(matFilePath);

if isa(matFile.N2Function,'function_handle') == true
    wvt = WVTransformHydrostatic([Lx Ly Lz],[Nx Ny Nz], N2=matFile.N2Function, latitude=latitude, rho0=rho0, ...
        dLnN2=dLnN2,PFinv=PFinv,QGinv=QGinv,PF=PF,QG=QG,h=h,P=P,Q=Q,z=z,shouldAntialias=shouldAntialias);
else
    wvt = WVTransformHydrostatic([Lx Ly Lz],[Nx Ny Nz], rho=matFile.rhoFunction, latitude=latitude, rho0=rho0, ...
        dLnN2=dLnN2,PFinv=PFinv,QGinv=QGinv,PF=PF,QG=QG,h=h,P=P,Q=Q,z=z,shouldAntialias=shouldAntialias);
end

wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);

end