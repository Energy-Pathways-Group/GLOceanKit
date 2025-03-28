function du = diffXIntoArray(self,u,du,options)
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
    self        WVFastTransformDoublyPeriodicFFTW
    u (:,:,:)   double
    du (:,:,:)   double
    options.n (1,1)     double = 1
end

% Method 0: Simplest, not fastest (2.21s-2.24s)
du = self.dftX.transformBack(self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftX.transformForward(u));

% Method 1: works fine (2.10s)
% self.dftXComplexBuffer = self.dftX.transformForwardIntoArray(u,self.dftXComplexBuffer);
% self.dftXComplexBuffer = self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftXComplexBuffer;
% du = self.dftX.transformBackIntoArray(self.dftXComplexBuffer,du);

% Method 2: works fine (1.85s-1.88s)
self.dftXComplexBuffer = self.dftX.transformForwardIntoArray(u,self.dftXComplexBuffer);
self.dftXComplexBuffer = self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftXComplexBuffer;
[self.dftXComplexBuffer,du] = self.dftX.transformBackIntoArrayDestructive(self.dftXComplexBuffer,du);

end