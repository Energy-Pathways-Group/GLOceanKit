function operations = operationForDynamicalVariable(self,variableName,options)
arguments(Input)
    self WVTransform {mustBeNonempty}
end
arguments (Input,Repeating)
    variableName char
end
arguments (Input)
    options.flowComponent WVFlowComponent = []
end
arguments (Output)
    operations WVOperation
end

% We want to do the minimum computation required, which means we need to
% eliminate tranforming zeros. We will not consider the optimization where
% one wave mask is zero and the other is not, since that is a relatively
% small gain.
if isempty(options.flowComponent)
    maskCase = 'all';
else
    flowComponent = options.flowComponent;
    if any(flowComponent.maskA0(:)) && ~any(flowComponent.maskAp(:) | flowComponent.maskAm(:))
        % geostrophic only
        maskCase = 'geostrophic';
    elseif any(flowComponent.maskA0(:)) && any(flowComponent.maskAp(:) | flowComponent.maskAm(:))
        % geostrophic, and (by exlusion) the wave parts are different
        maskCase = 'geostrophic-wave';
    elseif ~any(flowComponent.maskA0(:)) && any(flowComponent.maskAp(:) | flowComponent.maskAm(:))
        % no geostrophic, and but some wave component
        maskCase = 'wave';
    elseif ~any(flowComponent.maskA0(:) | flowComponent.maskAp(:) | flowComponent.maskAm(:))
        error('Your mask eliminates all flow!');
    else
        error('This case should not be reachable!');
    end
end

% as far as I can tell, there's no way to actually preallocating this in a
% way that Matlab finds helpful.
operations = WVOperation.empty(length(variableName),0);

for iOp = 1:length(variableName)
    name = variableName{iOp};
    switch name
        case 'u'
            switch maskCase
                case 'all'
                    varAnnotation = WVVariableAnnotation(name,{'x','y','z'},'m/s', 'x-component of the fluid velocity');
                    varAnnotation.attributes('standard_name') = 'eastward_sea_water_velocity';
                case {'geostrophic','wave','geostrophic-wave'}
                    varAnnotation = WVVariableAnnotation(strcat(name,'_',flowComponent.abbreviatedName),{'x','y','z'},'m/s', strcat('x-component of the fluid velocity, ',flowComponent.name,' component'));
            end
            switch maskCase
                case 'all'
                    f = @(wvt) wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt,A0=wvt.UA0.*wvt.A0t);
                case {'geostrophic'}
                    f = @(wvt) wvt.transformToSpatialDomainWithF(A0=flowComponent.maskA0.*wvt.UA0.*wvt.A0t);
                case 'wave'
                    f = @(wvt) wvt.transformToSpatialDomainWithF(Apm=flowComponent.maskAp.*wvt.UAp.*wvt.Apt + flowComponent.maskAm.*wvt.UAm.*wvt.Amt);
                case 'geostrophic-wave'
                    f = @(wvt) wvt.transformToSpatialDomainWithF(Apm=flowComponent.maskAp.*wvt.UAp.*wvt.Apt + flowComponent.maskAm.*wvt.UAm.*wvt.Amt,A0=flowComponent.maskA0.*wvt.UA0.*wvt.A0t);
            end
            
    end
    operations(end+1) = WVOperation(varAnnotation.name,varAnnotation,f);
end

end