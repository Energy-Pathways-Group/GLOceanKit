function Finv = FinvMatrix(wvt)
% transformation matrix $$F^{-1}$$
%
% A matrix that transforms a vector from vertical mode space to physical
% space.
%
% - Topic: Operations — Transformations
% - Declaration: Finv = FinvMatrix(wvt)
% - Returns Finv: A matrix with dimensions [Nz Nj]
arguments
    wvt         WVTransform
end

Finv = shiftdim(wvt.P,1) .* wvt.PFinv;

end