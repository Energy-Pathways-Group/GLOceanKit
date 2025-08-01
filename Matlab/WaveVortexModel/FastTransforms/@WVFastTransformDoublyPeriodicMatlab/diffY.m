function du = diffY(self,u,options)
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
    self        WVFastTransformDoublyPeriodicMatlab
    u (:,:,:)   double
    options.n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*shiftdim(self.wvg.l_dft,-1)).^options.n .*fft(u,self.wvg.Ny,2), self.wvg.Ny, 2,'symmetric');

end