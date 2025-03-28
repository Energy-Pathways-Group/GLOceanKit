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
    self        WVFastTransformDoublyPeriodicFFTW
    u (:,:,:)   double
    options.n (1,1)     double = 1
end

% Method 1: works fine
self.dftYComplexBuffer = self.dftY.transformForwardIntoArray(u,self.dftYComplexBuffer);
self.dftYComplexBuffer = self.dftY.scaleFactor * ((self.dy).^options.n) .* self.dftYComplexBuffer;
du = double(zeros(self.dftY.realSize));
[self.dftYComplexBuffer,du] = self.dftY.transformBackIntoArrayDestructive(self.dftYComplexBuffer,du);

% Method 2: now passes, but is slower
% [self.dftYComplexBuffer,self.dftRealBuffer] = self.dftY.transformBackIntoArrayDestructive(self.dftYComplexBuffer,self.dftRealBuffer);
% du = self.dftRealBuffer;
end