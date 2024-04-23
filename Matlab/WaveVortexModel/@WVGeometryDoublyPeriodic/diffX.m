function du = diffX(wvg,u,n)
% differentiate a spatial variable in the x-direction
%
% Performs spectral differentiation on variable u.
%
% - Topic: Operations — Differentiation
% - Declaration: du = diffX(u,n)
% - Parameter u: variable with dimensions $$(x,y,z)$$
% - Parameter n: (optional) order of differentiation $$\frac{d^n}{dx^n}$$ (default 1)
% - Returns du: differentiated variable in the spatial domain
arguments
    wvg         WVGeometryDoublyPeriodic
    u (:,:,:)   double
    n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*wvg.k).^n .* fft(u,wvg.Nx,1), wvg.Nx, 1,'symmetric');

end