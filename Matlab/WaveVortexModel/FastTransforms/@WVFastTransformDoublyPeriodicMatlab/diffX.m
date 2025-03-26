function du = diffX(self,u,options)
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
    self        WVFastTransformDoublyPeriodicMatlab
    u (:,:,:)   double
    options.n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*self.wvg.k_dft).^options.n .* fft(u,self.wvg.Nx,1), self.wvg.Nx, 1,'symmetric');

end