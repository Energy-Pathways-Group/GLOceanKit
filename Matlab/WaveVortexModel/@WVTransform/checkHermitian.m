function A = checkHermitian(A)
% Check if the matrix is Hermitian. Report errors.
%
% This function makes assumptions about the structure of the matrix.
%
% The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
% numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
% conjugates. Also, the negative k wavenumbers for l=0.
%
% - Topic: Utility function
% - Declaration: A = checkHermitian( A )
% - Returns A: matrix the same size as the input matrix
    M = size(A,1);
    N = size(A,2);
    K = size(A,3);

    for k=1:K
        for i=M:-1:1
            for j=N:-1:1
                ii = mod(M-i+1, M) + 1;
                jj = mod(N-j+1, N) + 1;
                if A(i,j,k) ~= conj(A(ii,jj,k))
                    fprintf('(i,j,k)=(%d,%d,%d) is not conjugate with (%d,%d,%d)\n',i,j,k,ii,jj,k)
                end
            end
        end
    end
end