function du = diffY(wvt,u,n)
arguments
    wvt         WaveVortexTransform
    u (:,:,:)   double
    n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*shiftdim(self.l,-1)).^n .*fft(u,wvt.Ny,2), wvt.Ny, 2,'symmetric');

end