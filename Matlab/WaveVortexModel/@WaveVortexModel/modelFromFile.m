function model = modelFromFile(path,options)
    arguments
        path char {mustBeFile}
        options.restartIndex (1,1) double {mustBePositive} = Inf
        options.shouldDoubleResolution double {mustBeMember(options.shouldDoubleResolution,[0 1])} = 0 
    end

    wvt = WaveVortexTransform.transformFromFile(path,iTime=options.restartIndex);

    ncfile = NetCDFFile(path);
    if isKey(ncfile.attributes,'NonlinearFluxOperation')
        nlFluxClassName = ncfile.attributes('NonlinearFluxOperation');
        nlFlux = feval(strcat(nlFluxClassName,'.nonlinearFluxFromFile'),ncfile,wvt);
    end
    if options.shouldDoubleResolution == 1
        wvt = wvt.transformWithDoubleResolution();
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