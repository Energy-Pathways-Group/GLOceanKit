%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate a 3D matrix to be Hermitian, except at k=l=0
function A = generateHermitianRandomMatrix( size, options )
arguments
    size (1,:) double
    options.shouldExcludeNyquist (1,1) double=1
    options.allowMeanPhase (1,1) double=0
end

nX = size(1); nY = size(2);
if length(size) > 2
    nZ = size(3);
else
    nZ = 1;
end
A = WaveVortexTransform.makeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
if options.shouldExcludeNyquist == 1
    mask = ~WaveVortexTransform.NyquistWavenumbers(A(:,:,1));
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