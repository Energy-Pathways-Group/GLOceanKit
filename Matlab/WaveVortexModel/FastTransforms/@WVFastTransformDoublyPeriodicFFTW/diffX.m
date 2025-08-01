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
    self        WVFastTransformDoublyPeriodicFFTW
    u (:,:,:)   double
    options.n (1,1)     double = 1
end

% Matlab builtin gets to 1.52s-1.55s

% Method 0: Simplest, not fastest (2.21s-2.24s)
% du = self.dftX.transformBack(self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftX.transformForward(u));

% Method 1: works fine (1.96s-2.00s)
% self.dftXComplexBuffer = self.dftX.transformForwardIntoArray(u,self.dftXComplexBuffer);
% self.dftXComplexBuffer = self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftXComplexBuffer;
% du = double(zeros(self.dftX.realSize));
% du = self.dftX.transformBackIntoArray(self.dftXComplexBuffer,du);

% Method 2: works fine (1.75s-1.8s)
self.dftXComplexBuffer = self.dftX.transformForwardIntoArray(u,self.dftXComplexBuffer);
self.dftXComplexBuffer = self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftXComplexBuffer;
du = double(zeros(self.dftX.realSize));
[self.dftXComplexBuffer,du] = self.dftX.transformBackIntoArrayDestructive(self.dftXComplexBuffer,du);

% Method 3: now passes, but is slower (1.85s-2.04s)
% self.dftXComplexBuffer = self.dftX.transformForwardIntoArray(u,self.dftXComplexBuffer);
% self.dftXComplexBuffer = self.dftX.scaleFactor * ((self.dx).^options.n) .* self.dftXComplexBuffer;
% [self.dftXComplexBuffer,self.dftRealBuffer] = self.dftX.transformBackIntoArrayDestructive(self.dftXComplexBuffer,self.dftRealBuffer);
% du = self.dftRealBuffer;
end