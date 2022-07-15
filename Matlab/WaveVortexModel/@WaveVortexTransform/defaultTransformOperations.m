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

transformOperations(end+1) = TransformOperation('t',StateVariable('t',{}, 's', 'time of observations'),@(wvt) wvt.t);

transformVar = StateVariable('A0',{'k','l','j'},'m', 'geostrophic coefficients at reference time t0');
transformVar.isComplex = 1;
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('A0', transformVar,@(wvt) wvt.Ap);

transformVar = StateVariable('Ap',{'k','l','j'},'m/s', 'positive wave coefficients at reference time t0');
transformVar.isComplex = 1;
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('Ap', transformVar,@(wvt) wvt.A0);

transformVar = StateVariable('Am',{'k','l','j'},'m/s', 'negative wave coefficients at reference time t0');
transformVar.isComplex = 1;
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('Am', transformVar,@(wvt) wvt.Am);

transformVar = StateVariable('internalWaveEnergyPlus',{},'m3/s2', 'total energy, internal waves, positive');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('internalWaveEnergyPlus', transformVar,@(wvt) wvt.internalWaveEnergyPlus);

transformVar = StateVariable('internalWaveEnergyMinus',{},'m3/s2', 'total energy, internal waves, minus');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('internalWaveEnergyMinus', transformVar,@(wvt) wvt.internalWaveEnergyMinus);

transformVar = StateVariable('waveEnergy',{},'m3/s2', 'total energy, waves');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('waveEnergy', transformVar,@(wvt) wvt.waveEnergy);

transformVar = StateVariable('inertialEnergyBaroclinic',{},'m3/s2', 'total energy, inertial oscillations, baroclinic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('inertialEnergyBaroclinic',transformVar,@(wvt) wvt.inertialEnergyBaroclinic);

transformVar = StateVariable('inertialEnergyBarotropic',{},'m3/s2', 'total energy, inertial oscillations, barotropic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('inertialEnergyBarotropic',transformVar,@(wvt) wvt.inertialEnergyBarotropic);

transformVar = StateVariable('inertialEnergy',{},'m3/s2', 'total energy, inertial oscillations');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('inertialEnergy', transformVar,@(wvt) wvt.inertialEnergy);

transformVar = StateVariable('geostrophicEnergyBaroclinic',{},'m3/s2', 'total energy, geostrophic, baroclinic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('geostrophicEnergyBaroclinic',transformVar,@(wvt) wvt.geostrophicEnergyBaroclinic);

transformVar = StateVariable('geostrophicEnergyBarotropic',{},'m3/s2', 'total energy, geostrophic, barotropic');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('geostrophicEnergyBarotropic',transformVar,@(wvt) wvt.geostrophicEnergyBarotropic);

transformVar = StateVariable('geostrophicEnergy',{},'m3/s2', 'total energy, inertial oscillations');
transformVar.isVariableWithLinearTimeStep = 0;
transformVar.isVariableWithNonlinearTimeStep = 1;
transformOperations(end+1) = TransformOperation('geostrophicEnergy', transformVar,@(wvt) wvt.geostrophicEnergy);

outputVar(1) = StateVariable('Apt',{'k','l','j'},'m/s', 'positive wave coefficients at time t');
outputVar(1).isComplex = 1;

outputVar(2) = StateVariable('Amt',{'k','l','j'},'m/s', 'negative wave coefficients at time t');
outputVar(2).isComplex = 1;

outputVar(3) = StateVariable('A0t',{'k','l','j'},'m', 'geostrophic coefficients at time t');
outputVar(3).isComplex = 1;

f = @(wvt) wvt.waveVortexCoefficientsAtTimeT();
transformOperations(end+1) = TransformOperation('ApAmA0',outputVar,f);

outputVar = StateVariable('u',{'x','y','z'},'m/s', 'x-component of the fluid velocity');
f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt + wvt.UA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('u',outputVar,f);

outputVar = StateVariable('v',{'x','y','z'},'m/s', 'y-component of the fluid velocity');
f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.VAp.*wvt.Apt + wvt.VAm.*wvt.Amt + wvt.VA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('v',outputVar,f);

outputVar = StateVariable('w',{'x','y','z'},'m/s', 'z-component of the fluid velocity');
f = @(wvt) wvt.transformToSpatialDomainWithG(wvt.WAp.*wvt.Apt + wvt.WAm.*wvt.Amt);
transformOperations(end+1) = TransformOperation('w',outputVar,f);

outputVar = StateVariable('p',{'x','y','z'},'kg/m/s2', 'pressure anomaly');
f = @(wvt) wvt.rho0*wvt.g*wvt.transformToSpatialDomainWithF(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('p',outputVar,f);

outputVar = StateVariable('eta',{'x','y','z'},'m', 'isopycnal deviation');
f = @(wvt) wvt.transformToSpatialDomainWithG(wvt.NAp.*wvt.Apt + self.NAm.*wvt.Amt + self.NA0.*wvt.A0t);
transformOperations(end+1) = TransformOperation('eta',outputVar,f);

outputVar = StateVariable('qgpv',{'x','y','z'},'1/s', 'quasigeostrophic potential vorticity');
f = @(wvt) -wvt.transformToSpatialDomainWithF( (wvt.Omega .* wvt.Omega / (wvt.h * wvt.f0)) .*wvt.A0t);
transformOperations(end+1) = TransformOperation('qgpv',outputVar,f);

fluxVar(1) = StateVariable('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
fluxVar(2) = StateVariable('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
fluxVar(3) = StateVariable('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');
transformOperations(end+1) = TransformOperation('nonlinearFlux',fluxVar,@(wvt) wvt.nonlinearFlux);

end