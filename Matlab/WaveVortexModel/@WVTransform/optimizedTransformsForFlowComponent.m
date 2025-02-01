function [transformToSpatialDomainWithF,transformToSpatialDomainWithG,mask,isMasked] = optimizedTransformsForFlowComponent(primaryFlowComponents,flowComponent)
% returns optimized transforms that avoid unnecessary computation
arguments (Input)
    primaryFlowComponents WVFlowComponent = WVFlowComponent.empty(0,0)
    flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
end

% These will be empty if +operationForKnownVariable is called when building
% documentation.
if isempty(flowComponent) && isempty(primaryFlowComponents)
    transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt),A0=A0(wvt));
    transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt),A0=A0(wvt));
    mask.Ap = 1;
    mask.Am = 1;
    mask.A0 = 1;
    isMasked = false;
    return
end

% What are the default masks. They won't always be 1, because QG models
% won't have Ap, Am
unmask.Ap = 0;
unmask.Am = 0;
unmask.A0 = 0;
for i = 1:length(primaryFlowComponents)
    unmask.Ap = unmask.Ap | any(primaryFlowComponents(i).maskAp(:));
    unmask.Am = unmask.Am | any(primaryFlowComponents(i).maskAm(:));
    unmask.A0 = unmask.A0 | any(primaryFlowComponents(i).maskA0(:));
end

% The goal here to to take the masks that we are given, which are of
% spectralMatrixSize, and reduce them to the scalar 1 or 0, if we can.
if isempty(flowComponent)
    mask = unmasked;
else
    % There can be zeros in a mask if there aren't any degrees of freedom
    % there, so we need to check for that.
    fullAp = zeros(self.spectralMatrixSize);
    fullAm = zeros(self.spectralMatrixSize);
    fullA0 = zeros(self.spectralMatrixSize);
    for i = 1:length(primaryFlowComponents)
        fullAp = fullAp + primaryFlowComponents(i).maskAp;
        fullAm = fullAm + primaryFlowComponents(i).maskAm;
        fullA0 = fullA0 + primaryFlowComponents(i).maskA0;
    end

    if ~any(flowComponent.maskAp(:))
        mask.Ap = 0;
    elseif isequal( flowComponent.maskAp & fullAp, fullAp )
        mask.Ap = 1;
    else
        mask.Ap = flowComponent.maskAp;
    end

    if ~any(flowComponent.maskAm(:))
        mask.Am = 0;
    elseif isequal( flowComponent.maskAm & fullAm, fullAm )
        mask.Am = 1;
    else
        mask.Am = flowComponent.maskAm;
    end

    if ~any(flowComponent.maskA0(:))
        mask.A0 = 0;
    elseif isequal( flowComponent.maskA0 & fullA0, fullA0 )
        mask.A0 = 1;
    else
        mask.A0 = flowComponent.maskA0;
    end
end

isMasked = ~(isscalar(mask.Ap) && isscalar(mask.Am) && isscalar(mask.A0) && mask.Ap == unmasked.Ap && mask.Am == unmasked.Am && mask.A0 == unmasked.A0);

% Various cases where we can eliminate significant computational time.
if (isscalar(mask.Ap) && isscalar(mask.Am) && mask.Ap == 0 && mask.Am == 0)
    if isscalar(mask.A0) && mask.A0 == 1
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(A0=A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(A0=A0(wvt));
    elseif isscalar(mask.A0) && mask.A0 == 0
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF();
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG();
    else
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(A0=mask.A0.*A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(A0=mask.A0.*A0(wvt));
    end
elseif (isscalar(mask.Ap) && isscalar(mask.Am) && mask.Ap == 1 && mask.Am == 1)
    if isscalar(mask.A0) && mask.A0 == 1
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt),A0=A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt),A0=A0(wvt));
    elseif isscalar(mask.A0) && mask.A0 == 0
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt));
    else
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=Ap(wvt)+Am(wvt),A0=mask.A0.*A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=Ap(wvt)+Am(wvt),A0=mask.A0.*A0(wvt));
    end
else
    if isscalar(mask.A0) && mask.A0 == 1
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=mask.Ap.*Ap(wvt)+mask.Am.*Am(wvt),A0=A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=mask.Ap.*Ap(wvt)+mask.Am.*Am(wvt),A0=A0(wvt));
    elseif isscalar(mask.A0) && mask.A0 == 0
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=mask.Ap.*Ap(wvt)+mask.Am.*Am(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=mask.Ap.*Ap(wvt)+mask.Am.*Am(wvt));
    else
        transformToSpatialDomainWithF = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithF(Apm=mask.Ap.*Ap(wvt)+mask.Am.*Am(wvt),A0=mask.A0.*A0(wvt));
        transformToSpatialDomainWithG = @(wvt,Ap,Am,A0) wvt.transformToSpatialDomainWithG(Apm=mask.Ap.*Ap(wvt)+mask.Am.*Am(wvt),A0=mask.A0.*A0(wvt));
    end
end

end