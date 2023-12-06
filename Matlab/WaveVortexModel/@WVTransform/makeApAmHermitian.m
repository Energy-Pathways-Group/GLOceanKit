function [Ap,Am] = makeApAmHermitian(self,Ap,Am)
% Forces the Ap/Am matrices to have the correct symmetries
%
% This function is NOT a true "Make Hermitian" function because it
% the Ap/Am matrices do not require k=l=0 to be real.
%
% If conjugateDimension == 2, then the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
% numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
% conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
% numbers are set to zero to avoid complications.
%
% - Topic: Utility function
% - Declaration: [Ap,Am] = makeApAmHermitian(Ap,Am)
% - Returns Ap: matrix the same size as the input matrix
% - Returns Am: matrix the same size as the input matrix

K = size(Ap,1);
L = size(Ap,2);
if self.conjugateDimension == 1
    % The order of the for-loop is chosen carefully here.
    for iK=1:(K/2+1)
        for iL=1:L
            if iK == 1 && iL > L/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                continue;
            else
                [Ap,Am] = setConjugate(Ap,Am,iK,iL,K,L);
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
                [Ap,Am] = setConjugate(Ap,Am,iK,iL,K,L);
            end
        end
    end
else
    error('invalid conjugate dimension')
end
end

function [Ap,Am] = setConjugate(Ap,Am,iK,iL,K,L)
icK = mod(K-iK+1, K) + 1;
icL = mod(L-iL+1, L) + 1;
if iK == icK && iL == icL % self-conjugate terms
    % This is not normally what you'd do for an FFT matrix, but
    % we're being Ap/Am are slightly different
    if iK == 1 && iL == 1
        Am(iK,iL,:) = conj(Ap(iK,iL,:));
    else
        Ap(iK,iL,:) = 0;
        Am(iK,iL,:) = 0;
    end
elseif iL == L/2+1 % Kill the Nyquist, because its never resolved
    Ap(iK,iL,:) = 0;
    Am(iK,iL,:) = 0;
else 
    Ap(icK,icL,:) = conj(Ap(iK,iL,:));
    Am(icK,icL,:) = conj(Am(iK,iL,:));
end
end