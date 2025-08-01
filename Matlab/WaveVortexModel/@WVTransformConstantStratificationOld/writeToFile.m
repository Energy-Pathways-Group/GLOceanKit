function [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)
% Output the `WVTransformConstantStratification` instance to file.
%
% - Topic: Write to file
% - Declaration: ncfile = writeToFile(netcdfFile,variables,options)
% - Parameter path: path to write file
% - Parameter variables: strings of variable names.
% - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.
% - Parameter shouldAddDefaultVariables: (optional) boolean indicating whether or not add default variables `A0`,`Ap`,`Am`,`t`. Default 1.
arguments
    wvt WVTransform {mustBeNonempty}
    path char {mustBeNonempty}
end
arguments (Repeating)
    variables char
end
arguments
    options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
    options.shouldAddDefaultVariables double {mustBeMember(options.shouldAddDefaultVariables,[0 1])} = 1
    options.shouldUseClassicNetCDF double {mustBeMember(options.shouldUseClassicNetCDF,[0 1])} = 1
end
variables = union(variables,{'N0'});
[ncfile,matFilePath] = writeToFile@WVTransform(wvt,path,variables{:},shouldAddDefaultVariables=options.shouldAddDefaultVariables,shouldOverwriteExisting=options.shouldOverwriteExisting,shouldUseClassicNetCDF=options.shouldUseClassicNetCDF);

end