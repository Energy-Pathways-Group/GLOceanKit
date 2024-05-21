function Ginv = GinvMatrix(wvt)
% transformation matrix $$G^{-1}$$
%
% A matrix that transforms a vector from vertical mode space to physical
% space.
%
% - Topic: Operations — Transformations
% - Declaration: Ginv = GinvMatrix(wvt)
% - Returns Finv: A matrix with dimensions [Nz Nj]
arguments
    wvt         WVTransform
end

Ginv = shiftdim(wvt.Q0,1) .* wvt.QG0inv;

end