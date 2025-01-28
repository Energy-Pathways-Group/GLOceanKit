function addPrimaryFlowComponent(self,primaryFlowComponent)
% add a primary flow component, automatically added to the flow
% components
%
% - Topic: Flow components
% - Declaration: addPrimaryFlowComponent(primaryFlowComponent)
% - Parameter primaryFlowComponent: one or more WVPrimaryFlowComponent objects
arguments
    self WVTransform {mustBeNonempty}
    primaryFlowComponent (1,:) WVPrimaryFlowComponent {mustBeNonempty}
end
for i=1:length(primaryFlowComponent)
    self.primaryFlowComponentNameMap{primaryFlowComponent(i).shortName} = primaryFlowComponent(i);
    self.flowComponentNameMap{primaryFlowComponent(i).shortName} = primaryFlowComponent(i);

    self.maskApPrimary = self.maskApPrimary | primaryFlowComponent(i).maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
    self.maskAmPrimary = self.maskAmPrimary | primaryFlowComponent(i).maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Am);
    self.maskA0Primary = self.maskA0Primary | primaryFlowComponent(i).maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0);

    self.maskApConj = self.maskApConj | primaryFlowComponent(i).maskOfConjugateModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
    self.maskAmConj = self.maskAmConj | primaryFlowComponent(i).maskOfConjugateModesForCoefficientMatrix(WVCoefficientMatrix.Am);
    self.maskA0Conj = self.maskA0Conj | primaryFlowComponent(i).maskOfConjugateModesForCoefficientMatrix(WVCoefficientMatrix.A0);
end
end