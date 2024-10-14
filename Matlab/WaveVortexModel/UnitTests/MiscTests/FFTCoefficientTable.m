wvg = WVGeometryDoublyPeriodic([2*pi 2*pi], [16 16]);


function Aklz = wvToDFT(wvg,Azkl)
Aklz = nan(wvg.Nx*wvg.Ny,1);
for iK=1:wvg.Nkl_wv
    Aklz(wvg.dftPrimaryIndices(iK)) = Azkl(iK);
end
Aklz = reshape(Aklz,[wvg.Nx wvg.Ny]);
end

k_dft = wvToDFT(wvg,wvg.k_wv);
l_dft = wvToDFT(wvg,wvg.l_wv);

k_dft = fftshift(fftshift(k_dft,1),2);
l_dft = fftshift(fftshift(l_dft,1),2);