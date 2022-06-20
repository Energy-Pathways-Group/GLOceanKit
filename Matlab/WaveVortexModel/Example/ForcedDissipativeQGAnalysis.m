wvt = WaveVortexTransform.transformFromFile('ForcedDissipativeQG-256.nc',iTime=25);

% m^2/s^3
EkT = wvt.transformToRadialWavenumber((wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)));
dk = (wvt.kRadial(2)-wvt.kRadial(1));

figure
plot(wvt.kRadial,EkT/dk), xlog, ylog, hold on
ylabel('m^3/s^2')
xlabel('1/m')
title('horizontal velocity spectrum, spectral density')
vlines([k_f,k_r],'g--')


[K,L,~] = ndgrid(wvt.k,wvt.l,wvt.j);
Kh = sqrt(K.*K + L.*L);
kRadial = wvt.kRadial;
AA = ~(wvt.MaskForAliasedModes(jFraction=1));
energyFlux = zeros(length(kRadial),1);
for iK=1:length(kRadial)
    A0Mask = AA;
    A0Mask(Kh > kRadial(iK)-dk/2 & Kh < kRadial(iK)+dk/2) = 0;
    [Ep,Em,E0] = wvt.energyFluxWithMasks(zeros(size(wvt.A0)),zeros(size(wvt.A0)),A0Mask);
    energyFlux(iK) = sum(E0(:));
end
figure, plot(kRadial,energyFlux/wvt.h), xlog
title('energy flux into a given wavenumber')
ylabel('m^2/s^3')