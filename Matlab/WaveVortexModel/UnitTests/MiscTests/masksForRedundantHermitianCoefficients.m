function [Ap,Am] = masksForRedundantHermitianCoefficients(self)
% Returns a matrix with 1s at the 'redundant' Hermitian indices.
%
% Returns a matrix with 1s at the 'redundant' Hermitian indices.
%
% - Topic: Utility function
% - Declaration: [Ap,Am] = masksForRedundantHermitianCoefficients(self)
% - Returns Ap: matrix of size [Nk, Nl, Nj] containing 1s and 0s
% - Returns Am: matrix of size [Nk, Nl, Nj] containing 1s and 0s

Ap = zeros(self.Nk,self.Nl,self.Nj);
Am = zeros(self.Nk,self.Nl,self.Nj);

K = size(Ap,1);
L = size(Ap,2);
if self.conjugateDimension == 1
    % The order of the for-loop is chosen carefully here.
    for iK=1:(K/2+1)
        for iL=1:L
            if iK == 1 && iL > L/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                continue;
            else
                [Ap,Am] = setConjugateToUnity(Ap,Am,iK,iL,K,L);
            end
        end
    end
elseif self.conjugateDimension == 2
    % The order of the for-loop is chosen carefully here.
    for iL=1:(L/2+1)
        for iK=1:K
            if iL == 1 && iK > K/2 % avoid letting l=0, k=Nx/2+1 terms set themselves again
                continue;
            else
                [Ap,Am] = setConjugateToUnity(Ap,Am,iK,iL,K,L);
            end
        end
    end
else
    error('invalid conjugate dimension')
end
end

function [Ap,Am] = setConjugateToUnity(Ap,Am,iK,iL,K,L)
icK = mod(K-iK+1, K) + 1;
icL = mod(L-iL+1, L) + 1;
if iK == icK && iL == icL % self-conjugate terms
    % This is not normally what you'd do for an FFT matrix, but
    % we're being Ap/Am are slightly different
    if iK == 1 && iL == 1
        Am(iK,iL,:) = 1;
    else
        Ap(iK,iL,:) = 0;
        Am(iK,iL,:) = 0;
    end
elseif iL == L/2+1 % Kill the Nyquist, because its never resolved
    Ap(iK,iL,:) = 0;
    Am(iK,iL,:) = 0;
else 
    Ap(icK,icL,:) = 1;
    Am(icK,icL,:) = 1;
end
end