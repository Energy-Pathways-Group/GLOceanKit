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

requiredVariables = {'x','y','z','Lx','Ly','Lz','j','latitude','t0','t','rho0','dLnN2','shouldAntialias'};
requiredVariables = union(requiredVariables,{'PF0inv','QG0inv','PF0','QG0','P0','Q0','h_0'});
requiredVariables = union(requiredVariables,{'PFpmInv','QGpmInv','PFpm','QGpm','Ppm','Qpm','h_pm','QGwg'});
if ~all(isKey(ncfile.variableWithName,requiredVariables) | isKey(ncfile.attributes,requiredVariables))
    error('This files is missing required variables or attributes to initialize a WVTransform.')
end

[x,y,z,Lx,Ly,Lz,rho0,latitude,dLnN2,shouldAntialias] = ncfile.readVariables('x','y','z','Lx','Ly','Lz','rho0','latitude','dLnN2','shouldAntialias');
Nx = length(x);
Ny = length(y);
Nz = length(z);
[PF0inv,QG0inv,PF0,QG0,P0,Q0,h_0] = ncfile.readVariables('PF0inv','QG0inv','PF0','QG0','P0','Q0','h_0');
[PFpmInv,QGpmInv,PFpm,QGpm,Ppm,Qpm,h_pm,QGwg] = ncfile.readVariables('PFpmInv','QGpmInv','PFpm','QGpm','Ppm','Qpm','h_pm','QGwg');

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
    wvt = WVTransformBoussinesq([Lx Ly Lz],[Nx Ny Nz], N2=matFile.N2Function, latitude=latitude, rho0=rho0, ...
        dLnN2=dLnN2,z=z,shouldAntialias=shouldAntialias, ...
        PF0inv=PF0inv,QG0inv=QG0inv,PF0=PF0,QG0=QG0,h_0=h_0,P0=P0,Q0=Q0, ...
        PFpmInv=PFpmInv,QGpmInv=QGpmInv,PFpm=PFpm,QGpm=QGpm,h_pm=h_pm,Ppm=Ppm,Qpm=Qpm,QGwg=QGwg);
else
    wvt = WVTransformBoussinesq([Lx Ly Lz],[Nx Ny Nz], rho=matFile.rhoFunction, latitude=latitude, rho0=rho0, ...
        dLnN2=dLnN2,z=z,shouldAntialias=shouldAntialias, ...
        PF0inv=PF0inv,QG0inv=QG0inv,PF0=PF0,QG0=QG0,h_0=h_0,P0=P0,Q0=Q0, ...
        PFpmInv=PFpmInv,QGpmInv=QGpmInv,PFpm=PFpm,QGpm=QGpm,h_pm=h_pm,Ppm=Ppm,Qpm=Qpm,QGwg=QGwg);
end

wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);

end