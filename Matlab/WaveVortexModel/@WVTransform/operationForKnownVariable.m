function operations = operationForKnownVariable(variableName,options)
% This is one of two functions that returns operations for computing
% standard variables.
arguments (Input,Repeating)
    variableName char
end
arguments (Input)
    options.spectralDimensionNames = {'j','kl'}
    options.spatialDimensionNames = {'x','y','z'}
    options.primaryFlowComponents WVFlowComponent = WVFlowComponent.empty(0,0)
    options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
end
arguments (Output)
    operations WVOperation
end

[transformToSpatialDomainWithF,transformToSpatialDomainWithG,mask,isMasked] = WVTransform.optimizedTransformsForFlowComponent(options.primaryFlowComponents,options.flowComponent);

operations = WVOperation.empty(length(variableName),0);

knownMaskableVariables = ["u","v","w","eta","p","psi","qgpv","energy"];
d = setdiff(variableName,knownMaskableVariables);
if isMasked && ~isempty(d)
    error("The variables %s cannot be masked.",strjoin(d,", "));
end

for iOp = 1:length(variableName)
    name = variableName{iOp};
    switch name
        case 'u'
            varAnnotation = WVVariableAnnotation(name,options.spatialDimensionNames,'m/s', 'x-component of the fluid velocity');
            varAnnotation.attributes('standard_name') = 'eastward_sea_water_velocity';
            f = @(wvt) transformToSpatialDomainWithF(wvt,@(wvt) wvt.UAp.*wvt.Apt,@(wvt) wvt.UAm.*wvt.Amt,@(wvt) wvt.UA0.*wvt.A0t);

        case 'v'
            varAnnotation = WVVariableAnnotation(name,options.spatialDimensionNames,'m/s', 'y-component of the fluid velocity');
            varAnnotation.attributes('standard_name') = 'northward_sea_water_velocity';
            f = @(wvt) transformToSpatialDomainWithF(wvt,@(wvt) wvt.VAp.*wvt.Apt,@(wvt) wvt.VAm.*wvt.Amt,@(wvt) wvt.VA0.*wvt.A0t);

        case 'w'
            varAnnotation = WVVariableAnnotation(name,options.spatialDimensionNames,'m/s', 'z-component of the fluid velocity');
            varAnnotation.attributes('standard_name') = 'upwardward_sea_water_velocity';
            f = @(wvt) transformToSpatialDomainWithG(wvt,@(wvt) wvt.WAp.*wvt.Apt,@(wvt) wvt.WAm.*wvt.Amt,@(wvt) 0);

        case 'eta'
            varAnnotation = WVVariableAnnotation(name,options.spatialDimensionNames,'m', 'approximate isopycnal deviation');
            f = @(wvt) transformToSpatialDomainWithG(wvt,@(wvt) wvt.NAp.*wvt.Apt,@(wvt) wvt.NAm.*wvt.Amt,@(wvt) wvt.NA0.*wvt.A0t);

        case 'p'
            varAnnotation =  WVVariableAnnotation(name,options.spatialDimensionNames,'kg/m/s2', 'pressure anomaly');
            f = @(wvt) wvt.rho0*wvt.g*transformToSpatialDomainWithF(wvt,@(wvt) wvt.NAp.*wvt.Apt,@(wvt) wvt.NAm.*wvt.Amt,@(wvt) wvt.PA0.*wvt.A0t);

        case 'psi'
            varAnnotation =  WVVariableAnnotation(name,options.spatialDimensionNames,'m^2/s', 'geostrophic streamfunction');
            varAnnotation.isVariableWithLinearTimeStep = 0;
            varAnnotation.isVariableWithNonlinearTimeStep = 1;
            f = @(wvt) (wvt.g/wvt.f) *transformToSpatialDomainWithF(wvt,@(wvt) 0,@(wvt) 0,@(wvt) wvt.PA0.*wvt.A0t);

        case 'qgpv'
            varAnnotation =  WVVariableAnnotation(name,options.spatialDimensionNames,'1/s', 'quasigeostrophic potential vorticity');
            varAnnotation.isVariableWithLinearTimeStep = 0;
            varAnnotation.isVariableWithNonlinearTimeStep = 1;
            f = @(wvt) transformToSpatialDomainWithF(wvt,@(wvt) 0,@(wvt) 0,@(wvt) wvt.A0_QGPV_factor .*wvt.A0t);

        case 'energy'
            varAnnotation = WVVariableAnnotation('energy',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
            varAnnotation.isVariableWithLinearTimeStep = 0;
            varAnnotation.isVariableWithNonlinearTimeStep = 1;
            f = @(wvt) sum( wvt.Apm_TE_factor(:).*( mask.Ap(:).*abs(wvt.Ap(:)).^2 + mask.Am(:).*abs(wvt.Am(:)).^2 ) + wvt.A0_TE_factor(:).*( mask.A0(:).*abs(wvt.A0(:)).^2) );

        case 'uvMax'
            varAnnotation = WVVariableAnnotation('uvMax',{},'m s^{-1}', 'max horizontal fluid speed');
            f = @(wvt) max(max(max( sqrt( (wvt.u).^2 + (wvt.v).^2 ) )));

        case 'wMax'
            varAnnotation = WVVariableAnnotation('wMax',{},'m s^{-1}', 'max vertical fluid speed');
            f = @(wvt) max(max(max( abs(wvt.w)  )));

        case 'rho_e'
            varAnnotation = WVVariableAnnotation('rho_e',options.spatialDimensionNames,'kg/m3', 'excess density');
            f = @(wvt) (wvt.rho0/wvt.g) * shiftdim(wvt.N2,-2) .* wvt.eta;

        case 'rho_total'
            varAnnotation = WVVariableAnnotation('rho_total',options.spatialDimensionNames,'kg/m3', 'total potential density');
            f = @(wvt) reshape(wvt.rho_nm,1,1,[]) + wvt.rho_e;

        case 'zeta_z'
            varAnnotation = WVVariableAnnotation('zeta_z',options.spatialDimensionNames,'1/s', 'vertical component of relative vorticity');
            varAnnotation.attributes('short_name') = 'ocean_relative_vorticity';
            f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);

        case 'ssu'
            varAnnotation = WVVariableAnnotation('ssu',{'x','y'},'m/s', 'x-component of the fluid velocity at the surface',detailedDescription='- topic: State Variables');
            f = @(wvt) wvt.u(:,:,end);

        case 'ssv'
            varAnnotation = WVVariableAnnotation('ssv',{'x','y'},'m/s', 'y-component of the fluid velocity at the surface',detailedDescription='- topic: State Variables');
            f = @(wvt) wvt.v(:,:,end);

        case 'ssh'
            varAnnotation = WVVariableAnnotation('ssh',{'x','y'},'m', 'sea-surface height');
            f = @(wvt) wvt.p(:,:,end)/(wvt.rho0*wvt.g);

        otherwise
            error('There is no variable named %s.',name)
    end

    if isMasked
        newName = strcat(varAnnotation.name,'_',options.flowComponent.abbreviatedName);
        newDescription = append(varAnnotation.description,', ',options.flowComponent.name,' component');
        varAnnotation = WVVariableAnnotation(newName,varAnnotation.dimensions,varAnnotation.units, newDescription);
    end

    operations(iOp) = WVOperation(varAnnotation.name,varAnnotation,f);
end

end