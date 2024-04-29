function [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = generateRandomFlowState(self)
% Random flow state, separated out by solution type.
%
% Generate a complete set of wave-vortex coefficients with variance at all
% physically realizable solution states.
%
% This is useful for testing that the transformation matrices are really
% complete and that the energy is diagonalizable.
% 
% Adding the solution types together, gives a complete state.
% Ap = ApIO + ApIGW;
% Am = AmIO + AmIGW;
% A0 = A0G + A0G0 + A0rhobar;
%
% - Declaration: [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = generateRandomFlowState()
% - Topic: Utility function

randCoeffs = @() randn(self.spectralMatrixSize) + sqrt(-1)*randn(self.spectralMatrixSize);

solutionGroup = WVInternalGravityWaveSolutionGroup(self);
ApIGW = randCoeffs() .* solutionGroup.maskForCoefficientMatrix(WVCoefficientMatrix.Ap);
AmIGW = randCoeffs() .* solutionGroup.maskForCoefficientMatrix(WVCoefficientMatrix.Am);

solutionGroup = WVGeostrophicSolutionGroup(self);
A0G = 6e-2*randCoeffs() .* solutionGroup.maskForCoefficientMatrix(WVCoefficientMatrix.A0);
A0G0 = zeros(self.spectralMatrixSize);
A0G0(self.J==0) = A0G(self.J==0);
A0G(self.J==0) = 0;

solutionGroup = WVInertialOscillationSolutionGroup(self);
ApIO = randCoeffs() .* solutionGroup.maskForCoefficientMatrix(WVCoefficientMatrix.Ap);
AmIO = conj(ApIO);

solutionGroup = WVMeanDensityAnomalySolutionGroup(self);
A0rhobar = randn(self.spectralMatrixSize) .* solutionGroup.maskForCoefficientMatrix(WVCoefficientMatrix.A0);
end