- Topic: Domain Attributes — Grid – Spectral

This is the usual definition of wavenumbers (or frequencies) for the Fourier transform,
```matlab
dl = 1/self.Ly;  
l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
```

The negative wavenumbers follow the positive wavenumbers. Matlab's built in function `fftshift` is useful for making these monotonically increasing.
