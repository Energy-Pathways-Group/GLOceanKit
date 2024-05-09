function self = addPrimaryFlowComponents(self)
flowComponent = WVGeostrophicComponent(self);
self.addPrimaryFlowComponent(flowComponent);
[self.A0Z,self.A0N] = flowComponent.geostrophicSpectralTransformCoefficients;
[self.UA0,self.VA0,self.NA0,self.PA0] = flowComponent.geostrophicSpatialTransformCoefficients;

flowComponent = WVMeanDensityAnomalyComponent(self);
self.addPrimaryFlowComponent(flowComponent);
A0Nmda = flowComponent.meanDensityAnomalySpectralTransformCoefficients;
NA0mda = flowComponent.meanDensityAnomalySpatialTransformCoefficients;
self.A0N = self.A0N + A0Nmda;
self.NA0 = self.NA0 + NA0mda;
self.PA0 = self.PA0 + NA0mda;

flowComponent = WVInternalGravityWaveComponent(self);
self.addPrimaryFlowComponent(flowComponent);
[self.ApmD,self.ApmN] = flowComponent.internalGravityWaveSpectralTransformCoefficients;
[self.UAp,self.VAp,self.WAp,self.NAp] = flowComponent.internalGravityWaveSpatialTransformCoefficients;

flowComponent = WVInertialOscillationComponent(self);
self.addPrimaryFlowComponent(flowComponent);
[UAp_io,VAp_io] = flowComponent.inertialOscillationSpatialTransformCoefficients;
self.UAp = self.UAp + UAp_io;
self.VAp = self.VAp + VAp_io;

self.UAm = conj(self.UAp);
self.VAm = conj(self.VAp);
self.WAm = self.WAp;
self.NAm = -self.NAp;

self.iOmega = sqrt(-1)*self.Omega;

self.Apm_TE_factor = zeros(self.spectralMatrixSize);
self.A0_TE_factor = zeros(self.spectralMatrixSize);
self.A0_QGPV_factor = zeros(self.spectralMatrixSize);
self.A0_TZ_factor = zeros(self.spectralMatrixSize);
for name = keys(self.primaryFlowComponentNameMap)
    flowComponent = self.primaryFlowComponentNameMap(name{1});
    self.Apm_TE_factor = self.Apm_TE_factor + flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.Ap);
    self.A0_TE_factor = self.A0_TE_factor + flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.A0);
    self.A0_QGPV_factor = self.A0_QGPV_factor + flowComponent.qgpvFactorForA0;
    self.A0_TZ_factor = self.A0_TZ_factor + flowComponent.enstrophyFactorForA0;
end
end