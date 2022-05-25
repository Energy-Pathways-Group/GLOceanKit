function du = diffX(wvt,u,n)
arguments
    wvt         WaveVortexTransform
    u (:,:,:)   double
    n (1,1)     double = 1
end

du = ifft( (sqrt(-1)*self.k).^n .* fft(u,wvt.Nx,1), wvt.Nx, 1,'symmetric');

end