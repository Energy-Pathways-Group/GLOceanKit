function wvt = waveVortexTransformFromFile(path,options)
% Initialize a WVTransform instance from an existing file
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

if isKey(ncfile.attributes,'WVTransform')
    wvtClassName = ncfile.attributes('WVTransform');
    wvt = feval(strcat(wvtClassName,'.waveVortexTransformFromFile'),path,'iTime',options.iTime);
end

if isKey(ncfile.attributes,'WVNonlinearFluxOperation')
    nlFluxClassName = ncfile.attributes('WVNonlinearFluxOperation');
    nlFlux = feval(strcat(nlFluxClassName,'.nonlinearFluxFromFile'),ncfile,wvt);
    wvt.nonlinearFluxOperation = nlFlux;
end

end