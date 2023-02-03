function wvt = waveVortexTransformFromFile(path,options)
% Initialize a WVTransformConstantStratification instance from an existing file
%
% - Topic: Initialization
% - Declaration: wvt = waveVortexTransformFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter iTime: (optional) time index to initialize from (default 1)
arguments
    path char {mustBeFile}
    options.iTime (1,1) double {mustBePositive} = 1
end

ncfile = NetCDFFile(path);

requiredVariables = {'x','y','z','Lx','Ly','Lz','latitude','N0','rho0'};
if ~all(isKey(ncfile.variableWithName,requiredVariables) | isKey(ncfile.attributes,requiredVariables))
    error('This files is missing required variables or attributes to initialize a WVTransform.')
end

[x,y,z,Lx,Ly,Lz,latitude,N0,rho0] = ncfile.readVariables('x','y','z','Lx','Ly','Lz','latitude','N0','rho0');
Nx = length(x);
Ny = length(y);
Nz = length(z);

wvt = WVTransformConstantStratification([Lx Ly Lz],[Nx Ny Nz],N0=N0,latitude=latitude,rho0=rho0);

wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);

end