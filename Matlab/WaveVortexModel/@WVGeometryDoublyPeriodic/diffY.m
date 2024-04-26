function du = diffY(wvg,u,n)
% differentiate a spatial variable in the y-direction
%
% Performs spectral differentiation on variable u.
%
% - Topic: Operations — Differentiation
% - Declaration: du = diffY(u,n)
% - Parameter u: variable with dimensions $$(x,y,z)$$
% - Parameter n: (optional) order of differentiation $$\frac{d^n}{dy^n}$$ (default 1)
% - Returns du: differentiated variable in the spatial domain
arguments
    wvg         WVGeometryDoublyPeriodic
    u (:,:,:)   double
    n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*shiftdim(wvg.l_dft,-1)).^n .*fft(u,wvg.Ny,2), wvg.Ny, 2,'symmetric');

end