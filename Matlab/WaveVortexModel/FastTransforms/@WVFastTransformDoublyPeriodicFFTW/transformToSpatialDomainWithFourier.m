function u = transformToSpatialDomainWithFourier(self,u_bar)
arguments
    self WVFastTransformDoublyPeriodic
    u_bar 
end

% Matlab builtin gets to 2.65s-2.70s

% Method 1: 2.63s-2.66s
self.dftXYComplexBuffer(self.wvg.dftPrimaryIndex) = u_bar;
self.dftXYComplexBuffer(self.wvg.dftConjugateIndex) = conj(u_bar(self.wvg.wvConjugateIndex));
u = self.dftXY.transformBack(self.dftXYComplexBuffer);

% Method 2: works (same timing as method 1)
% self.dftXYComplexBuffer(self.wvg.dftPrimaryIndex) = u_bar;
% self.dftXYComplexBuffer(self.wvg.dftConjugateIndex) = conj(u_bar(self.wvg.wvConjugateIndex));
% u = double(zeros(self.dftX.realSize));
% u = self.dftXY.transformBackIntoArray(self.dftXYComplexBuffer,u);

% Method 3: (3.43-3.45s)
% doesn't get copied, so gets stomped on later.
% self.dftXYComplexBuffer(self.wvg.dftPrimaryIndex) = u_bar;
% self.dftXYComplexBuffer(self.wvg.dftConjugateIndex) = conj(u_bar(self.wvg.wvConjugateIndex));
% self.dftRealBuffer = self.dftXY.transformBackIntoArray(self.dftXYComplexBuffer,self.dftRealBuffer);
% u = self.dftRealBuffer;

% Method 4: (2.46s-2.51s) [Doesn't pass unit tests... but why?]
% self.dftXYComplexBuffer(self.wvg.dftPrimaryIndex) = u_bar;
% self.dftXYComplexBuffer(self.wvg.dftConjugateIndex) = conj(u_bar(self.wvg.wvConjugateIndex));
% u = double(zeros(self.dftXY.realSize));
% [self.dftXYComplexBuffer,u] = self.dftXY.transformBackIntoArrayDestructive(self.dftXYComplexBuffer,u);


% [self.dftXYComplexBuffer,self.dftRealBuffer] = self.dftXY.transformBackIntoArrayDestructive(self.dftXYComplexBuffer,self.dftRealBuffer);
% u = self.dftRealBuffer;

% self.complexBuffer(self.wvg.dftPrimaryIndex) = u_bar;
% self.complexBuffer(self.wvg.dftConjugateIndex) = conj(u_bar(self.wvg.wvConjugateIndex));
% u = ifft(ifft(self.dftXYComplexBuffer,self.wvg.Nx,1),self.wvg.Ny,2,'symmetric')*(self.wvg.Nx*self.wvg.Ny);
end