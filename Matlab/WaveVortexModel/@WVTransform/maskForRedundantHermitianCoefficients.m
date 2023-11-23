function A = maskForRedundantHermitianCoefficients(self)
% Returns a matrix with 1s at the 'redundant' hermiation indices.
%
% This function makes assumptions about the structure of the matrix.
%
% Returns a matrix with 1s at the 'redundant' hermiation indices.
%
% - Topic: Masks
% - Declaration: A = maskForRedundantHermitianCoefficients(self)
% - Returns A: matrix of size [Nk, Nl, Nj] containing 
A = zeros(self.Nk,self.Nl,self.Nj);

M = self.Nk;
N = self.Nl;
K = self.Nj;

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