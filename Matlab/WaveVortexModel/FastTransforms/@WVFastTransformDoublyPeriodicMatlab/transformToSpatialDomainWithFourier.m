function u = transformToSpatialDomainWithFourier(self,u_bar)
self.complexBuffer(self.wvg.dftPrimaryIndex) = u_bar;
self.complexBuffer(self.wvg.dftConjugateIndex) = conj(u_bar(self.wvg.wvConjugateIndex));
u = ifft(ifft(self.complexBuffer,self.wvg.Nx,1),self.wvg.Ny,2,'symmetric')*(self.wvg.Nx*self.wvg.Ny);
end