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

    wvt = WaveVortexTransform.waveVortexTransformFromFile(path,iTime=options.restartIndex);

    ncfile = NetCDFFile(path);
    if isKey(ncfile.attributes,'WVNonlinearFluxOperation')
        nlFluxClassName = ncfile.attributes('WVNonlinearFluxOperation');
        nlFlux = feval(strcat(nlFluxClassName,'.nonlinearFluxFromFile'),ncfile,wvt);
    end
    if options.shouldDoubleResolution == 1
        wvt = wvt.waveVortexTransformWithDoubleResolution();
        if exist('nlFlux','var')
            nlFlux = nlFlux.nonlinearFluxWithDoubleResolution(wvt);
        end
    end

    if exist('nlFlux','var')
        model = WaveVortexModel(wvt,nonlinearFlux=nlFlux);
    else
        model = WaveVortexModel(wvt);
    end
    
% if there's existing model output, use that output interval
%     time = ncread(existingModelOutput,'t');
%     if length(time)>1
%         self.outputInterval = time(2)-time(1);
%     end
end