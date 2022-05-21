%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Takes a Hermitian matrix and resets it back to the real amplitude.
function [A,phi,linearIndex] = ExtractNonzeroWaveProperties(Matrix)
M = size(Matrix,1);
N = size(Matrix,2);
K = size(Matrix,3);

A = [];
phi = [];
linearIndex = [];

% The order of the for-loop is chosen carefully here.
for k=1:K
    for j=1:(N/2+1)
        for i=1:M
            ii = mod(M-i+1, M) + 1;
            jj = mod(N-j+1, N) + 1;
            waveAmp = 0; wavePhase = 0;
            if i == ii && j == jj
                % self-conjugate term
                if i == 1 && j == 1
                    waveAmp = abs(Matrix(i,j,k));
                    wavePhase = angle(Matrix(i,j,k));
                else
                    continue;
                end
            elseif j == N/2+1 % Kill the Nyquist, rather than fix it.
                waveAmp = abs(Matrix(i,j,k));
                wavePhase = angle(Matrix(i,j,k));
            else % we are letting l=0, k=Nx/2+1 terms set themselves again, but that's okay 
%                 A(ii,jj,k) = conj(A(i,j,k));
                if j == 1 && i > M/2
                    continue;
                end
                waveAmp = 2*abs(Matrix(i,j,k));
                wavePhase = angle(Matrix(i,j,k));
            end
            if waveAmp > 0
               A = cat(1,A,waveAmp);
               phi = cat(1,phi,wavePhase);
               linearIndex = cat(1,linearIndex,sub2ind(size(Matrix),i,j,k));
            end
        end
    end
end