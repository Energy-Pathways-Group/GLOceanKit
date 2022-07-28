function A = nyquistWavenumbers(A)
% Returns a matrix with 1s at the Nyquist frequencies.
%
% This function makes assumptions about the structure of the matrix.
%
% Returns a matrix of the same size as A with 1s at the Nyquist
% frequencies.
%
% - Topic: Utility function
% - Declaration: A = nyquistWavenumbers( A )
% - Returns A: matrix the same size as the input matrix
    A = zeros(size(A));
    Nx = size(A,1);
    Ny = size(A,2);
    if Nx > 1
        A(ceil(Nx/2)+1,:) = 1;
    end
    if Ny > 1
        A(:,ceil(Ny/2)+1) = 1;
    end
end