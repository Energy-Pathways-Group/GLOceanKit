function A = redundantHermitianCoefficients(A)
% Returns a matrix with 1s at the 'redundant' hermiation indices.
%
% This function makes assumptions about the structure of the matrix.
%
% Returns a matrix the same size as A with 1s at the 'redundant'
% hermiation indices.
%
% - Topic: Utility function
% - Declaration: A = redundantHermitianCoefficients( A )
% - Returns A: matrix the same size as the input matrix
    A = zeros(size(A));

    M = size(A,1);
    N = size(A,2);
    K = size(A,3);

    % The order of the for-loop is chosen carefully here.
    for k=1:K
        for j=1:(N/2+1)
            for i=1:M
                ii = mod(M-i+1, M) + 1;
                jj = mod(N-j+1, N) + 1;
                if i == ii && j == jj
                    if i == 1 && j == 1
                        continue;
                    else
                        A(i,j,k) = 0;
                    end
                elseif j == N/2+1
                    A(i,j,k) = 0;
                else
                    if j == 1 && i > M/2
                        continue;
                    end
                    A(ii,jj,k) = 1;
                end
            end
        end
    end
end