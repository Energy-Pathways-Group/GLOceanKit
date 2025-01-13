function [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)
% Output the `WVTransform` to file.
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
variables = union(variables,{'rho_nm','N2'});
[ncfile,matFilePath] = writeToFile@WVTransform(wvt,path,variables{:},shouldAddDefaultVariables=options.shouldAddDefaultVariables,shouldOverwriteExisting=options.shouldOverwriteExisting,shouldUseClassicNetCDF=options.shouldUseClassicNetCDF);

% To keep things tidy, lets put the transform in a separate group
group = ncfile.addGroup("WVTransformHydrostatic");
variables = {'dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h'};
for iVar=1:length(variables)
    if isKey(wvt.propertyAnnotationNameMap,variables{iVar})
        varAnnotation = wvt.propertyAnnotationNameMap(variables{iVar});
    elseif isKey(wvt.variableAnnotationWithName,variables{iVar})
        varAnnotation = wvt.variableAnnotationWithName(variables{iVar});
    end
    varAnnotation.attributes('units') = varAnnotation.units;
    varAnnotation.attributes('long_name') = varAnnotation.description;
    group.addVariable(varAnnotation.name,varAnnotation.dimensions,wvt.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
end

rhoFunction = wvt.rhoFunction;
N2Function = wvt.N2Function;
dLnN2Function = wvt.dLnN2Function;
date_created = ncfile.attributes('date_created');
save(matFilePath,'rhoFunction','N2Function','dLnN2Function','date_created');
fprintf('In addition to the NetCDF file, a .mat sidecar file was created at the same path.\n');
end
