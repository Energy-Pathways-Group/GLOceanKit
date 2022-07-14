function du = diffY(wvt,u,n)
% differentiate a spatial variable in the y-direction
%
% Performs spectral differentiation on variable u.
%
% - Topic: Operations — differentiation
% - Declaration: du = diffY(u,n)
% - Parameter u: variable with dimensions $$(x,y,z)$$
% - Parameter n: (optional) order of differentiation d^n/dy^n (default 1)
% - Returns du: differentiated variable in the spatial domain
arguments
    wvt         WaveVortexTransform
    u (:,:,:)   double
    n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*shiftdim(wvt.l,-1)).^n .*fft(u,wvt.Ny,2), wvt.Ny, 2,'symmetric');

end