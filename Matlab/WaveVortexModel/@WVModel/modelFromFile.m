function model = modelFromFile(path,options)
% Initialize a model from an existing file
%
% - Topic: Initialization
% - Declaration: model = modelFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter restartIndex: (optional) time index to initialize from (default 1)
% - Parameter shouldDoubleResolution: (optional) whether or not to double the resolution
    arguments
        path char {mustBeFile}
        options.restartIndex (1,1) double {mustBePositive} = Inf
        options.shouldDoubleResolution double {mustBeMember(options.shouldDoubleResolution,[0 1])} = 0 
    end

    wvt = WVTransform.waveVortexTransformFromFile(path,iTime=options.restartIndex);

    if options.shouldDoubleResolution == 1
        wvt = wvt.waveVortexTransformWithDoubleResolution();
    end

    ncfile = NetCDFFile(path);
    if ncfile.attributes('shouldUseLinearDynamics') == 0
        model = WVModel(wvt,nonlinearFlux=wvt.nonlinearFluxOperation);
    else
        model = WVModel(wvt);
    end
    
% if there's existing model output, use that output interval
%     time = ncread(existingModelOutput,'t');
%     if length(time)>1
%         self.outputInterval = time(2)-time(1);
%     end
end