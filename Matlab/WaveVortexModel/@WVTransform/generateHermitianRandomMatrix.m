function A = generateHermitianRandomMatrix( size, options )
% Generate a 3D matrix to be Hermitian, except at k=l=0
%
% This function makes assumptions about the structure of the matrix.
%
% Normally values at index (1,1,:) must be self-conjugate (and therefore
% real), but we added a flag to allow a phase, which is helpful for
% randomizing the phase of the inertial oscillations.
%
% - Topic: Utility function
% - Declaration: A = generateHermitianRandomMatrix( size, options )
% - Parameter size: size of the matrix to generate
% - Parameter shouldExcludeNyquist: optional (default 1) will set Nyquist frequencies to zero if true.
% - Parameter allowMeanPhase: optional (default 0) will all a compex component to values at index (1,1,:) if set to true
% - Returns A: Hermitian conjugate matrix of given size
arguments
    size (1,:) double
    options.shouldExcludeNyquist (1,1) double=1
    options.allowMeanPhase (1,1) double=0
end

nX = size(1); nY = size(2);
A = WVTransform.makeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
if options.shouldExcludeNyquist == 1
    mask = ~WVTransform.nyquistWavenumbers(A(:,:,1));
    A = mask.*A;
else
    A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
    A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
    A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
end

if options.allowMeanPhase == 1
    A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
else
    A(1,1,:) = 2*real(A(1,1,:)); % Double the zero frequency
end

end