function du = diffX(wvt,u,n)
% differentiate a spatial variable in the x-direction
%
% Performs spectral differentiation on variable u.
%
% - Topic: Operations — differentiation
% - Declaration: du = diffX(u,n)
% - Parameter u: variable with dimensions $$(x,y,z)$$
% - Parameter n: (optional) order of differentiation d^n/dx^n (default 1)
% - Returns du: differentiated variable in the spatial domain
arguments
    wvt         WaveVortexTransform
    u (:,:,:)   double
    n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*wvt.k).^n .* fft(u,wvt.Nx,1), wvt.Nx, 1,'symmetric');

end