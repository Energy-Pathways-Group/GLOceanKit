function F = FMatrix(wvt)
% transformation matrix $$F$$
%
% A matrix that transforms a vector from physical
% space to vertical mode space.
%
% - Topic: Operations — Transformations
% - Declaration: F = FMatrix(wvt)
% - Returns Finv: A matrix with dimensions [Nz Nj]
arguments
    wvt         WVTransform
end

F = wvt.PF ./ shiftdim(wvt.P,2);

end