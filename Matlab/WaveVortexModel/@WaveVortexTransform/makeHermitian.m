%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forces a 3D matrix to be Hermitian, ready for an FFT (internal)
%
% The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
% numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
% conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
% numbers are set to zero to avoid complications.
%
% This function is NOT a true "Make Hermitian" function because it
% doesn't force the k=l=0 to be real.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = makeHermitian(A)
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
                    % A(i,j,k) = real(A(i,j,k)); % self-conjugate term
                    % This is not normally what you'd do, but we're being
                    % tricky by later adding the conjugate
                    if i == 1 && j == 1
                        continue;
                    else
                        A(i,j,k) = 0;
                    end
                elseif j == N/2+1 % Kill the Nyquist, rather than fix it.
                    A(i,j,k) = 0;
                else % we are letting l=0, k=Nx/2+1 terms set themselves again, but that's okay
                    A(ii,jj,k) = conj(A(i,j,k));
                end
            end
        end
    end
end