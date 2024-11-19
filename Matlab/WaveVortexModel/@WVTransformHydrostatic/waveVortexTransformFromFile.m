function [wvt,ncfile] = waveVortexTransformFromFile(path,options)
% Initialize a WVTransformHydrostatic instance from an existing file
%
% This static method is called by WVTransform.waveVortexTransformFromFile
% and should not need to be called directly.
%
% - Topic: Initialization (Static)
% - Declaration: wvt = waveVortexTransformFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter iTime: (optional) time index to initialize from (default 1)
arguments (Input)
    path char {mustBeFile}
    options.iTime (1,1) double {mustBePositive} = 1
end
arguments (Output)
    wvt WVTransform
    ncfile NetCDFFile
end

ncfile = NetCDFFile(path);

requiredVariables = {'x','y','z','Lx','Ly','Lz','latitude','t0','rho0','dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h','shouldAntialias'};
if ~all(ncfile.hasVariableWithName(requiredVariables{:}))
    error('This files is missing required variables or attributes to initialize a WVTransform.')
end

[x,y,z,Lx,Ly,Lz,rho0,latitude,dLnN2,PF0inv,QG0inv,PF0,QG0,P0,Q0,h,shouldAntialias] = ncfile.readVariables('x','y','z','Lx','Ly','Lz','rho0','latitude','dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h','shouldAntialias');
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
        dLnN2=dLnN2,PFinv=PF0inv,QGinv=QG0inv,PF=PF0,QG=QG0,h=h,P=P0,Q=Q0,z=z,shouldAntialias=shouldAntialias);
else
    wvt = WVTransformHydrostatic([Lx Ly Lz],[Nx Ny Nz], rho=matFile.rhoFunction, latitude=latitude, rho0=rho0, ...
        dLnN2=dLnN2,PFinv=PF0inv,QGinv=QG0inv,PF=PF0,QG=QG0,h=h,P=P0,Q=Q0,z=z,shouldAntialias=shouldAntialias);
end

wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);

end