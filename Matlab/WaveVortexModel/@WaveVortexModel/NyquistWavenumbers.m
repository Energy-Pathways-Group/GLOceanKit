%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns a matrix the same size as A with 1s at the Nyquist
% frequencies.
function A = NyquistWavenumbers(A)
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