function [Fp,Fm,F0] = nonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)
% returns the flux of each coefficient as determined by the nonlinear flux operation
%
% This computes the nonlinear flux that results from a subset of flow
% constituents. The masks are applied to the coefficients Ap,Am,A0 before
% computing the nonlinear flux, $$\vec{u} \cdot \nabla \vec{u}$$. This
% function calls -nonlinearFluxWithGradientMasks.
%
% - Topic: Nonlinear flux and energy transfers
% - Declaration: [Fp,Fm,F0] = nonlinearFluxForFlowConstituents(Uconstituent,gradUconstituent)
% - Parameter Uconstituent: `WVFlowConstituent` type for $$\vec{u} \cdot \nabla \vec{u}$$
% - Parameter gradUconstituent: `WVFlowConstituent` type for $$\vec{u} \cdot \nabla \vec{u}$$
% - Returns Fp: flux into the Ap coefficients
% - Returns Fm: flux into the Am coefficients
% - Returns F0: flux into the A0 coefficients
arguments
    self WVTransform {mustBeNonempty}
    uFlowComponent WVFlowComponent
    gradUFlowComponent WVFlowComponent
end

[Fp,Fm,F0] = self.nonlinearFluxWithGradientMasks(uFlowComponent.maskAp,uFlowComponent.maskAm,uFlowComponent.maskA0,gradUFlowComponent.maskAp,gradUFlowComponent.maskAm,gradUFlowComponent.maskA0);
end