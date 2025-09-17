function [wvt,ncfile] = waveVortexTransformFromFile(path,options)
% Initialize a WVTransform instance from an existing file
%
% A WVTransform instance can be recreated from a NetCDF file and .mat
% sidecar file if the default variables were save to file. For example,
%
% ```matlab
% wvt = WVTransform.waveVortexTransformFromFile('cyprus-eddy.nc',iTime=Inf);
% ```
%
% will create a WVTransform instance populated with values from the last
% time point of the file. Note that this is a static function---it is a
% function defined on the class, not an instance variable---so requires we
% prepend `WVTransform.` The result of this function call is an instance
% variable.
%
% If you intend to read more than one time point from the save file, hold
% onto the NetCDFFile instance that is returned, and then call
% [`initFromNetCDFFile`](/classes/wvtransform/initfromnetcdffile.html). This
% avoids the relatively expensive operation recreating the WVTransform, and
% simply read the appropriate data from file.
%
% See also the users guide for [reading and writing to
% file](/users-guide/reading-and-writing-to-file.html).
%
% - Topic: Initialization
% - Declaration: [wvt,ncfile] = WVTransform.waveVortexTransformFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter iTime: (optional) time index to initialize from (default 1).
% - Returns wvt: an instance of a WVTransform subclass
% - Returns ncfile: a NetCDFFile instance pointing to the file
arguments (Input)
    path char {mustBeFile}
    options.iTime (1,1) double {mustBePositive} = 1
    options.shouldReadOnly logical = false
end
arguments (Output)
    wvt WVTransform
    ncfile NetCDFFile
end
ncfile = NetCDFFile(path,shouldReadOnly=options.shouldReadOnly);

if isKey(ncfile.attributes,'WVTransform')
    wvtClassName = ncfile.attributes('WVTransform');
    wvt = feval(strcat(wvtClassName,'.waveVortexTransformFromFile'),path,'iTime',options.iTime,'shouldReadOnly',options.shouldReadOnly);
else
    error("Unable to find the attribute WVTransform in this file, suggesting this was not created by the wave vortex model.");
end

% totalForcingGroups = ncfile.attributes('TotalForcingGroups');
% for iForce=1:totalForcingGroups
%     forceGroup = ncfile.groupWithName("forcing-"+iForce);
%     forcingClassName = forceGroup.attributes('WVForcing');
%     force = feval(strcat(forcingClassName,'.forcingFromFile'),forceGroup,wvt);
%     self.addForcing(force);
% end

end