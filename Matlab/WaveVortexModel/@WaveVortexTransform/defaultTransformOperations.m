function transformOperations = defaultTransformOperations()
% return array of TransformOperation instances initialized by default
%
% This function creates a number of standard StateVariables with associated
% TransformOperations.
%
% - Topic: Internal
% - Declaration: transformOperations = defaultTransformOperations()
% - Returns transformOperations: array of TransformOperation instances
transformOperations = TransformOperation.empty(0,0);

transformOperations(end+1) = TransformOperation('t',WVVariableAnnotation('t',{}, 's', 'time of observations'),@(wvt) wvt.t);

transformVar = WVVariableAnnotation('A0',{'k','l','j'},'m', 'geostrophic coefficients at reference time t0');
transformVar.isComplex = 1;
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('A0', transformVar,@(wvt) wvt.Ap);

transformVar = WVVariableAnnotation('Ap',{'k','l','j'},'m/s', 'positive wave coefficients at reference time t0');
transformVar.isComplex = 1;
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('Ap', transformVar,@(wvt) wvt.A0);

transformVar = WVVariableAnnotation('Am',{'k','l','j'},'m/s', 'negative wave coefficients at reference time t0');
transformVar.isComplex = 1;
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('Am', transformVar,@(wvt) wvt.Am);

transformVar = WVVariableAnnotation('totalEnergy',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('totalEnergy', transformVar,@(wvt) wvt.totalEnergy);

transformVar = WVVariableAnnotation('totalEnergySpatiallyIntegrated',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed in the spatial domain');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('totalEnergySpatiallyIntegrated', transformVar,@(wvt) wvt.totalEnergySpatiallyIntegrated);

transformVar = WVVariableAnnotation('totalHydrostaticEnergy',{},'m3/s2', 'horizontally-averaged depth-integrated energy *without w* computed in the spatial domain');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('totalHydrostaticEnergy', transformVar,@(wvt) wvt.totalHydrostaticEnergy);

transformVar = WVVariableAnnotation('internalWaveEnergyPlus',{},'m3/s2', 'total energy, internal waves, positive');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('internalWaveEnergyPlus', transformVar,@(wvt) wvt.internalWaveEnergyPlus);

transformVar = WVVariableAnnotation('internalWaveEnergyMinus',{},'m3/s2', 'total energy, internal waves, minus');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('internalWaveEnergyMinus', transformVar,@(wvt) wvt.internalWaveEnergyMinus);

transformVar = WVVariableAnnotation('waveEnergy',{},'m3/s2', 'total energy, waves');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('waveEnergy', transformVar,@(wvt) wvt.waveEnergy);

transformVar = WVVariableAnnotation('inertialEnergyBaroclinic',{},'m3/s2', 'total energy, inertial oscillations, baroclinic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('inertialEnergyBaroclinic',transformVar,@(wvt) wvt.inertialEnergyBaroclinic);

transformVar = WVVariableAnnotation('inertialEnergyBarotropic',{},'m3/s2', 'total energy, inertial oscillations, barotropic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('inertialEnergyBarotropic',transformVar,@(wvt) wvt.inertialEnergyBarotropic);

transformVar = WVVariableAnnotation('inertialEnergy',{},'m3/s2', 'total energy, inertial oscillations');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('inertialEnergy', transformVar,@(wvt) wvt.inertialEnergy);

transformVar = WVVariableAnnotation('geostrophicEnergyBaroclinic',{},'m3/s2', 'total energy, geostrophic, baroclinic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('geostrophicEnergyBaroclinic',transformVar,@(wvt) wvt.geostrophicEnergyBaroclinic);

transformVar = WVVariableAnnotation('geostrophicEnergyBarotropic',{},'m3/s2', 'total energy, geostrophic, barotropic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('geostrophicEnergyBarotropic',transformVar,@(wvt) wvt.geostrophicEnergyBarotropic);

transformVar = WVVariableAnnotation('geostrophicEnergy',{},'m3/s2', 'total energy, inertial oscillations');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('geostrophicEnergy', transformVar,@(wvt) wvt.geostrophicEnergy);

outputVar(1) = WVVariableAnnotation('Apt',{'k','l','j'},'m/s', 'positive wave coefficients at time t');
outputVar(1).isComplex = 1;

outputVar(2) = WVVariableAnnotation('Amt',{'k','l','j'},'m/s', 'negative wave coefficients at time t');
outputVar(2).isComplex = 1;

outputVar(3) = WVVariableAnnotation('A0t',{'k','l','j'},'m', 'geostrophic coefficients at time t');
outputVar(3).isComplex = 1;

f = @(wvt) wvt.waveVortexCoefficientsAtTimeT();
transformOperations(end+1) = TransformOperation('ApAmA0',outputVar,f);

outputVar = WVVariableAnnotation('uMax',{},'m s^{-1}', 'max horizontal fluid speed');
f = @(wvt) max(max(max( sqrt( (wvt.u).^2 + (wvt.v).^2 ) )));
transformOperations(end+1) = TransformOperation('uMax',outputVar,f);

outputVar = WVVariableAnnotation('u',{'x','y','z'},'m/s', 'x-component of the fluid velocity');
f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt + wvt.UA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('u',outputVar,f);

outputVar = WVVariableAnnotation('v',{'x','y','z'},'m/s', 'y-component of the fluid velocity');
f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.VAp.*wvt.Apt + wvt.VAm.*wvt.Amt + wvt.VA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('v',outputVar,f);

outputVar = WVVariableAnnotation('w',{'x','y','z'},'m/s', 'z-component of the fluid velocity');
f = @(wvt) wvt.transformToSpatialDomainWithG(wvt.WAp.*wvt.Apt + wvt.WAm.*wvt.Amt);
transformOperations(end+1) = TransformOperation('w',outputVar,f);

outputVar = WVVariableAnnotation('p',{'x','y','z'},'kg/m/s2', 'pressure anomaly');
f = @(wvt) wvt.rho0*wvt.g*wvt.transformToSpatialDomainWithF(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('p',outputVar,f);

outputVar = WVVariableAnnotation('eta',{'x','y','z'},'m', 'isopycnal deviation');
f = @(wvt) wvt.transformToSpatialDomainWithG(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('eta',outputVar,f);

outputVar = WVVariableAnnotation('qgpv',{'x','y','z'},'1/s', 'quasigeostrophic potential vorticity');
f = @(wvt) -wvt.transformToSpatialDomainWithF( (wvt.Omega .* wvt.Omega / (wvt.h * wvt.f)) .*wvt.A0t);
transformOperations(end+1) = TransformOperation('qgpv',outputVar,f);

fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap',detailedDescription='- topic: State Variables');
fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am',detailedDescription='- topic: State Variables');
fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0',detailedDescription='- topic: State Variables');
transformOperations(end+1) = TransformOperation('nonlinearFlux',fluxVar,@(wvt) wvt.nonlinearFlux);

end