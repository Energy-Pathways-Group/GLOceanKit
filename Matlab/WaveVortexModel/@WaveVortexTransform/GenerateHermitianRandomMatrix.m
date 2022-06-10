%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate a 3D matrix to be Hermitian, except at k=l=0
function A = GenerateHermitianRandomMatrix( size, shouldExcludeNyquist )
    nX = size(1); nY = size(2);
    if length(size) > 2
        nZ = size(3);
    else
        nZ = 1;
    end
    A = InternalWaveModel.makeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
    if shouldExcludeNyquist == 1
        mask = ~InternalWaveModel.NyquistWavenumbers(A(:,:,1));
        A = mask.*A;
    else
        A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
        A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
        A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
        A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
    end
end