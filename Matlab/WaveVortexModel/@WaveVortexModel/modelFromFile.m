function model = modelFromFile(path,options)
    arguments
        path char {mustBeFile}
        options.restartIndex (1,1) double {mustBePositive} = Inf
        options.shouldDoubleResolution double {mustBeMember(options.shouldDoubleResolution,[0 1])} = 0 
    end

    wvt = WaveVortexTransform.transformFromFile(path,iTime=options.restartIndex);
    if options.shouldDoubleResolution == 1
        wvt = wvt.transformWithDoubleResolution();
    end

    model = WaveVortexModel(wvt);
    
% if there's existing model output, use that output interval
%     time = ncread(existingModelOutput,'t');
%     if length(time)>1
%         self.outputInterval = time(2)-time(1);
%     end
end