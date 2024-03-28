function G = GMatrix(wvt)
% transformation matrix $$G$$
%
% A matrix that transforms a vector from physical
% space to vertical mode space.
%
% - Topic: Operations — Transformations
% - Declaration: G = GMatrix(wvt)
% - Returns Ginv: A matrix with dimensions [Nz Nj]
arguments
    wvt         WVTransform
end

G = wvt.QG ./ shiftdim(wvt.Q,2);

end