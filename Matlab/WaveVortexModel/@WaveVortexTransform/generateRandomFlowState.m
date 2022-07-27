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
ApIGW = WaveVortexTransform.generateHermitianRandomMatrix( size(self.Ap), shouldExcludeNyquist=1, allowMeanPhase=1 );
AmIGW = WaveVortexTransform.generateHermitianRandomMatrix( size(self.Ap), shouldExcludeNyquist=1, allowMeanPhase=1 );
A0G = 6e-2*WaveVortexTransform.generateHermitianRandomMatrix( size(self.Ap), shouldExcludeNyquist=1 );

ApIO = zeros(size(self.Ap));
AmIO = zeros(size(self.Ap));
A0G0 = zeros(size(self.Ap));
A0rhobar = zeros(size(self.Ap));

% inertial oscillations only exist at k=l=0
ApIO(1,1,:) = ApIGW(1,1,:);
AmIO(1,1,:) = conj(ApIGW(1,1,:));

% zero out all j=0, and k=l=0 values.
ApIGW(:,:,1) = 0;
ApIGW(1,1,:) = 0;
AmIGW(:,:,1) = 0;
AmIGW(1,1,:) = 0;

% barotropic geostrophic at all k and l>0, j=0
A0G0(:,:,1) = 0.1*A0G(:,:,1);
A0G0(1,1,1) = 0;

% mean density anomaly
A0rhobar(1,1,2:end) = real(A0G(1,1,2:end));

% zero out all j=0, and k=l=0 values.
A0G(1,1,:) = 0;
A0G(:,:,1) = 0;
end