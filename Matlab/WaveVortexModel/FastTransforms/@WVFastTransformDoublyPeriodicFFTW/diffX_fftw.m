function du = diffX_fftw(wvg,u,n)
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
    wvg         WVGeometryDoublyPeriodic
    u (:,:,:)   double
    n (1,1)     double = 1
end

wvg.fftw_complex_cache = wvg.dft.scaleFactor*(sqrt(-1)*wvg.k_hc).^n .* wvg.dft.transformForwardIntoArray(u,wvg.fftw_complex_cache);
% [wvg.fftw_complex_cache, wvg.fftw_real_cache] = wvg.dft.transformBackIntoArrayDestructive(wvg.fftw_complex_cache,wvg.fftw_real_cache);
wvg.fftw_real_cache = wvg.dft.transformBackIntoArrayDestructive(wvg.fftw_complex_cache,wvg.fftw_real_cache);
% [wvg.fftw_real_cache,wvg.fftw_complex_cache] = wvg.dft.transformBackIntoArrayDestructive(wvg.fftw_complex_cache,wvg.fftw_real_cache);

du = wvg.fftw_real_cache;

end