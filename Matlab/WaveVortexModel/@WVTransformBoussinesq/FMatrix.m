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

F = wvt.PF0 ./ shiftdim(wvt.P0,2);

end