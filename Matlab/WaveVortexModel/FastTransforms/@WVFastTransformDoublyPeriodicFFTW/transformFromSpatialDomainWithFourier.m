function u_bar = transformFromSpatialDomainWithFourier(self,u)
% transform from $$(x,y,z)$$ to $$(z,kl)$$ on the WV grid
%
% Performs a Fourier transform in the x and y direction. The
% resulting matrix is on the WV gride
%
% - Topic: Operations — Fourier transformation
% - Declaration: u_bar = transformFromSpatialDomainWithFourier(u)
% - Parameter u: a real-valued matrix of size [Nx Ny Nz]
% - Returns u_bar: a complex-valued matrix of size [Nz Nkl]

% self.dftXYComplexBuffer = self.dftXY.scaleFactor*self.dftXY.transformForwardIntoArray(u,self.dftXYComplexBuffer);

% Matlab builtin gets to 2.23s-2.25s

% Method 1: (2.7s)
% from performance testing it is clear that the multiplication step is
% causing a memory copy.
u_bar = self.dftXY.transformForward(u);
u_bar = u_bar*self.dftXY.scaleFactor;

% Method 2: (3.23s-3.30s)
% u_bar = complex(double(zeros(self.dftXY.complexSize)));
% u_bar = self.dftXY.scaleFactor*self.dftXY.transformForwardIntoArray(u,u_bar);


u_bar = reshape(u_bar(self.wvg.dftPrimaryIndex),[self.Nz self.wvg.Nkl]);

end