function operations = operationForDynamicalVariable(self,variableName,options)
% This function is designed with the following goals:
%   1) the definition of each variable is made *once*, with the intention
%   of minimizing mistakes and making it easier to read. The challenge with
%   the goals is then that masks have to be applied generically.
%   2) Masking an entire component (Apm or A0) to zero should avoid any
%   extra computation. This is important because we still want the QG flow
%   to retain its speed.
%
% The solution I came up with was that each variables must specify the
% recipe for computing itself as three function handles that tell it how to
% construct the three components Ap,Am,A0. The function then only applies
% the function handle when necessary.
arguments(Input)
    self WVTransform {mustBeNonempty}
end
arguments (Input,Repeating)
    variableName char
end
arguments (Input)
    options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
end
arguments (Output)
    operations WVOperation
end

% The goal here to to take the masks that we are given, which are of
% spectralMatrixSize, and reduce them to the scalar 1 or 0, if we can.
if isempty(options.flowComponent)
    % if the user didn't specify a mask, then all the masks are just 1.
    maskAp = 1;
    maskAm = 1;
    maskA0 = 1;
else

    % There can be zeros in a mask if there aren't any degrees of freedom
    % there, so we need to check for that.
    fullAp = zeros(self.spectralMatrixSize);
    fullAm = zeros(self.spectralMatrixSize);
    fullA0 = zeros(self.spectralMatrixSize);
    for name = keys(self.primaryFlowComponentNameMap)
        primaryFowComponent = self.primaryFlowComponentNameMap(name{1});
        fullAp = fullAp + primaryFowComponent.maskAp;
        fullAm = fullAm + primaryFowComponent.maskAm;
        fullA0 = fullA0 + primaryFowComponent.maskA0;
    end

    flowComponent = options.flowComponent;
    if ~any(flowComponent.maskAp(:))
        maskAp = 0;
    elseif isequal( flowComponent.maskAp & fullAp, fullAp )
        maskAp = 1;
    else
        maskAp = flowComponent.maskAp;
    end

    if ~any(flowComponent.maskAm(:))
        maskAm = 0;
    elseif isequal( flowComponent.maskAm & fullAm, fullAm )
        maskAm = 1;
    else
        maskAm = flowComponent.maskAm;
    end

    if ~any(flowComponent.maskA0(:))
        maskA0 = 0;
    elseif isequal( flowComponent.maskA0 & fullA0, fullA0 )
        maskA0 = 1;
    else
        maskA0 = flowComponent.maskA0;
    end
end
% true if at least one component is masked
isMasked = ~(isscalar(maskAp) && isscalar(maskAm) && isscalar(maskA0) && maskAp == 1 && maskAm == 1 && maskA0 == 1);

% Various cases where we can eliminate significant computational time.
if (isscalar(maskAp) && isscalar(maskAm) && maskAp == 0 && maskAm == 0)
    if isscalar(maskA0) && maskA0 == 1
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(A0=A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(A0=A0(wvt));
    elseif isscalar(maskA0) && maskA0 == 0
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF();
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG();
    else
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(A0=maskA0.*A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(A0=maskA0.*A0(wvt));
    end
elseif (isscalar(maskAp) && isscalar(maskAm) && maskAp == 1 && maskAm == 1)
    if isscalar(maskA0) && maskA0 == 1
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt),A0=A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt),A0=A0(wvt));
    elseif isscalar(maskA0) && maskA0 == 0
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt));
    else
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt),A0=maskA0.*A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt),A0=maskA0.*A0(wvt));
    end
else
    if isscalar(maskA0) && maskA0 == 1
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=maskAp.*Ap(wvt)+maskAm.*Am(wvt),A0=A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=maskAp.*Ap(wvt)+maskAm.*Am(wvt),A0=A0(wvt));
    elseif isscalar(maskA0) && maskA0 == 0
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=maskAp.*Ap(wvt)+maskAm.*Am(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=maskAp.*Ap(wvt)+maskAm.*Am(wvt));
    else
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=maskAp.*Ap(wvt)+maskAm.*Am,A0=maskA0.*A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=maskAp.*Ap(wvt)+maskAm.*Am,A0=maskA0.*A0(wvt));
    end
end

% as far as I can tell, there's no way to actually preallocating this in a
% way that Matlab finds helpful.
operations = WVOperation.empty(length(variableName),0);

for iOp = 1:length(variableName)
    name = variableName{iOp};
    switch name
        case 'u'
            varAnnotation = WVVariableAnnotation(name,{'x','y','z'},'m/s', 'x-component of the fluid velocity');
            varAnnotation.attributes('standard_name') = 'eastward_sea_water_velocity';
            f = @(wvt) transformToSpatialDomainWithF(wvt,@(wvt) wvt.UAp.*wvt.Apt,@(wvt) wvt.UAm.*wvt.Amt,@(wvt) wvt.UA0.*wvt.A0t);

        case 'v'
            varAnnotation = WVVariableAnnotation(name,{'x','y','z'},'m/s', 'y-component of the fluid velocity');
            varAnnotation.attributes('standard_name') = 'northward_sea_water_velocity';
            f = @(wvt) transformToSpatialDomainWithF(wvt,@(wvt) wvt.VAp.*wvt.Apt,@(wvt) wvt.VAm.*wvt.Amt,@(wvt) wvt.VA0.*wvt.A0t);

        case 'w'
            varAnnotation = WVVariableAnnotation(name,{'x','y','z'},'m/s', 'z-component of the fluid velocity');
            varAnnotation.attributes('standard_name') = 'upwardward_sea_water_velocity';
            f = @(wvt) transformToSpatialDomainWithG(wvt,@(wvt) wvt.WAp.*wvt.Apt,@(wvt) wvt.WAm.*wvt.Amt,@(wvt) 0);

        case 'eta'
            varAnnotation = WVVariableAnnotation(name,{'x','y','z'},'m', 'approximate isopycnal deviation');
            f = @(wvt) transformToSpatialDomainWithG(wvt,@(wvt) wvt.NAp.*wvt.Apt,@(wvt) wvt.NAm.*wvt.Amt,@(wvt) wvt.NA0.*wvt.A0t);

        case 'p'
            varAnnotation =  WVVariableAnnotation(name,{'x','y','z'},'kg/m/s2', 'pressure anomaly');
            f = @(wvt) wvt.rho0*wvt.g*transformToSpatialDomainWithF(wvt,@(wvt) wvt.NAp.*wvt.Apt,@(wvt) wvt.NAm.*wvt.Amt,@(wvt) wvt.PA0.*wvt.A0t);

        case 'psi'
            varAnnotation =  WVVariableAnnotation(name,{'x','y','z'},'m^2/s', 'geostrophic streamfunction');
            f = @(wvt) (wvt.g/wvt.f) *transformToSpatialDomainWithF(wvt,@(wvt) 0,@(wvt) 0,@(wvt) wvt.PA0.*wvt.A0t);

        case 'qgpv'
            varAnnotation =  WVVariableAnnotation(name,{'x','y','z'},'1/s', 'quasigeostrophic potential vorticity');
            f = @(wvt) transformToSpatialDomainWithF(wvt,@(wvt) 0,@(wvt) 0,@(wvt) wvt.A0_QGPV_factor .*wvt.A0t);

        case 'energy'
            varAnnotation = WVVariableAnnotation('energy',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
            f = @(wvt) sum( wvt.Apm_TE_factor(:).*( maskAp(:).*abs(wvt.Ap(:)).^2 + maskAm(:).*abs(wvt.Am(:)).^2 ) + wvt.A0_TE_factor(:).*( maskA0(:).*abs(wvt.A0(:)).^2) );

        otherwise
            error('There is no dynamical variable named %s.',name)
    end

    if isMasked == 1
        newName = strcat(varAnnotation.name,'_',options.flowComponent.abbreviatedName);
        newDescription = append(varAnnotation.description,', ',options.flowComponent.name,' component');
        varAnnotation = WVVariableAnnotation(newName,varAnnotation.dimensions,varAnnotation.units, newDescription);
    end

    operations(iOp) = WVOperation(varAnnotation.name,varAnnotation,f);
end

end