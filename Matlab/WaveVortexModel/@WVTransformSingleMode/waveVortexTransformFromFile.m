function wvt = waveVortexTransformFromFile(path,options)
% Initialize a WVTransformSingleMode instance from an existing file
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

requiredVariables = {'x','y','Lx','Ly','Lz','latitude'};
if ~all(isKey(ncfile.variableWithName,requiredVariables) | isKey(ncfile.attributes,requiredVariables))
    error('This files is missing required variables or attributes to initialize a WVTransform.')
end

[x,y,Lx,Ly,Lz,latitude] = ncfile.readVariables('x','y','Lx','Ly','Lz','latitude');
Nx = length(x);
Ny = length(y);

wvt = WVTransformSingleMode([Lx Ly],[Nx Ny],h=Lz,latitude=latitude);

wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);

end