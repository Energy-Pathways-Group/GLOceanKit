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
u_bar = fft(fft(u,self.wvg.Nx,1),self.wvg.Ny,2)/(self.wvg.Nx*self.wvg.Ny);
u_bar = reshape(u_bar(self.wvg.dftPrimaryIndex),[self.Nz self.wvg.Nkl]);
end