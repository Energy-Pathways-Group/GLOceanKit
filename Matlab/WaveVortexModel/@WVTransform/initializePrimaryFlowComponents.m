function self = initializePrimaryFlowComponents(self)

flowComponent = WVInternalGravityWaveComponent(self);
self.addPrimaryFlowComponent(flowComponent);
[self.ApmD,self.ApmN] = flowComponent.internalGravityWaveSpectralTransformCoefficients;
[self.UAp,self.VAp,self.WAp,self.NAp] = flowComponent.internalGravityWaveSpatialTransformCoefficients;

self.UAm = conj(self.UAp);
self.VAm = conj(self.VAp);
self.WAm = self.WAp;
self.NAm = -self.NAp;

self.iOmega = sqrt(-1)*self.Omega;

self.Apm_TE_factor = zeros(self.spectralMatrixSize);
self.A0_TE_factor = zeros(self.spectralMatrixSize);
self.A0_QGPV_factor = zeros(self.spectralMatrixSize);
self.A0_TZ_factor = zeros(self.spectralMatrixSize);
for name =self.primaryFlowComponentNames
    flowComponent = self.primaryFlowComponent(name{1});
    self.Apm_TE_factor = self.Apm_TE_factor + flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.Ap);
    self.A0_TE_factor = self.A0_TE_factor + flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.A0);
    self.A0_QGPV_factor = self.A0_QGPV_factor + flowComponent.qgpvFactorForA0;
    self.A0_TZ_factor = self.A0_TZ_factor + flowComponent.enstrophyFactorForA0;
end
end