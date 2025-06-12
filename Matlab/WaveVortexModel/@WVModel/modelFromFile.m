function model = modelFromFile(path)
% Initialize a model from an existing file
%
% A WVModel will be initialized from the specified path. The model will
% have this file designated as its outputFile. Integrating the model will
% thus write to this file.
%
% - Topic: Initialization
% - Declaration: model = modelFromFile(path)
% - Parameter path: path to a NetCDF file
    arguments
        path char {mustBeFile}
    end

    wvt = WVTransform.waveVortexTransformFromFile(path,iTime=Inf);
    model = WVModel(wvt);
    outputFile = WVModelOutputFile.modelOutputFileFromFile(NetCDFFile(path),model);
    model.addOutputFile(outputFile);
end